"""
    WAVSource(io)

A readble audio stream that reads from a WAV file represented by the given
io stream. The stream will accumulate metadata chunks as it reads through
the file, which you can access with the `metadata` function. Note that if the
file has any metadata chunks after the data chunk in the file, they won't be
available until the whole file is read.
"""
mutable struct WAVSource{IO, T} <: SampleSource
    io::IO
    format::WAVFormat
    opt::Dict{String, Vector{Vector{UInt8}}}
    file_bytesleft::Int
    data_bytesleft::Int
end

# This method reads the WAV header and everything up to the "data" chunk,
# leaving the stream ready to be read from.
function WAVSource(io)
    bytesleft = read_header(io)
    fmt, opt, bytesleft = find_format(io, bytesleft)
    eltype = type_from_format(fmt)
    # read into the stream until we find the data
    data_bytesleft, bytesleft = find_data(io, opt, bytesleft)

    # now we have a WAVSource ready to provide audio data
    WAVSource{typeof(io), eltype}(io, fmt, opt, bytesleft, data_bytesleft)
end

SampledSignals.nchannels(src::WAVSource) = Int(src.format.nchannels)
SampledSignals.samplerate(src::WAVSource) = float(src.format.sample_rate)
Base.eltype(::WAVSource{IO, T}) where {IO, T} = T

SampledSignals.metadata(src::WAVSource, key) = src.opt[key][1]
SampledSignals.metadata(src::WAVSource, key, ::Colon) = src.opt[key]
SampledSignals.metadata(src::WAVSource, key, idx) = src.opt[key][idx]

# These are the MATLAB compatible signatures
wavread(filename::AbstractString, fmt::AbstractString) = wavread(filename, format=fmt)
wavread(filename::AbstractString, n) = wavread(filename, subrange=n)
wavread(filename::AbstractString, n, fmt) = wavread(filename, subrange=n, format=fmt)

# implement SampledSignals's read method to plug into the stream conversion
# infrastructure
function SampledSignals.unsafe_read!(src::WAVSource, buf::Array,
                                     frameoffset, framecount)
    framebytes = src.format.block_align
    samplebytes = framebytes ÷ nchannels(src)
    bytestoread = if isnull(src.data_bytesleft)
            framecount * framebytes
        else
            min(framecount * framebytes, src.data_bytesleft)
        end

    # TODO: pre-allocate this buffer
    rawdata = read(src.io, bytestoread)
    bytesread = length(rawdata)
    if bytesread < bytestoread
        # got EOF while reading
        if(!isnull(src.file_bytesleft))
            warn("reached EOF but RIFF header indicated more data")
        end
        if(!isnull(src.data_bytesleft))
            warn("reached EOF but data chunk header indicated more data")
        end
        src.file_bytesleft -= 0
        src.data_bytesleft -= 0
    else
        src.file_bytesleft -= bytesread
        src.data_bytesleft -= bytesread
    end
    framesread = bytesread ÷ framebytes

    decode_buffer(buf, rawdata, src.format, frameoffset, framesread, samplebytes)

    if src.data_bytesleft < framebytes
        # we've reached the end of the data chunk
        if src.data_bytesleft > 0
            warn("\"data\" chunk had $(src.data_bytesleft) orphaned bytes")
            skip(src.io, src.data_bytesleft)
            src.file_bytesleft -= src.data_bytesleft
            src.data_bytesleft = 0
        end

        # grab any of the chunks that follow
        parse_tail(src.io, src.opt, src.file_bytesleft)
    end

    framesread
end

"""
    function decode_buffer(dest, src,
                           format, frameoffset, framecount,
                           samplebytes)
Copy the samples from `src` (the raw UInt8 buffer from the stream) into `dest`
(the array used for unsafe_read)
"""
function decode_buffer(dest::Array{T, N}, src,
                       format, frameoffset, framecount,
                       samplebytes) where {T, N}
    rawidx = 1
    for frame in (1:framecount) + frameoffset
        for ch in 1:size(dest, 2)
            dest[frame, ch] = decodewavbytes(T, src, rawidx, samplebytes)
            rawidx += samplebytes
        end
    end
end

# we need special handling for reading into 16-bit because we need to
# check whether we actually have an 8-bit companded stream
function decode_buffer(dest::Array{Fixed{Int16, 15}, N}, src,
                       format, frameoffset, framecount,
                       samplebytes) where N
    rawidx = 1
    if isformat(format, WAVE_FORMAT_PCM)
        # OK, just a normal 16-bit PCM
        for frame in (1:framecount) + frameoffset, ch in 1:size(dest, 2)
            dest[frame, ch] = decodewavbytes(Fixed{Int16, 15}, src, rawidx, 2)
            rawidx += 2
        end
    elseif isformat(format, WAVE_FORMAT_MULAW)
        for frame in (1:framecount) + frameoffset, ch in 1:size(dest, 2)
            dest[frame, ch] = decodemulaw(src[rawidx])
            rawidx += 1
        end
    elseif isformat(format, WAVE_FORMAT_ALAW)
        for frame in (1:framecount) + frameoffset, ch in 1:size(dest, 2)
            dest[frame, ch] = decodealaw(src[rawidx])
            rawidx += 1
        end
    else
        throw_fmt_error(format)
    end

    framecount
end


# decodings:
# 8-bit PCM: UInt8 -> PCM8Sample (needs offset - `reinterpret(PCM8Sample, x-0x80)`)
# 8-bit alaw: UInt8 -> PCM16Sample
# 8-bit mulaw: UInt8 -> PCM16Sample
# UInt16 -> PCM16Sample (reinterpret)
# UInt32 -> PCM32Sample (reinterpret)
# UInt64 -> PCM64Sample (reinterpret)
# UInt32 -> Float32 (reinterpret)
# UInt64 -> Float64 (reinterpret)

"""
    decodewavbytes(eltype, rawdata, rawidx, nbytes)

Decode raw bytes read from a wav file into host-endian samples of type `eltype`.
The data is pulled from the `rawdata` vector, starting at `rawidx`. the `nbytes`
argument specifies how many bytes the sample takes up in the buffer, and is only
used for PCM types that can be varying in size (e.g. 24-bit audio takes 3 bytes,
but the target type is Int32).

These functions should be type-stable and fast, so they can be called in an
inner loop.
"""
function decodewavbytes end

decodewavbytes(::Type{UInt8}, rawdata, rawidx, nbytes) = rawdata[rawidx]

# create a native unsigned integer type from the raw bytes
function decodewavbytes(T::Type{<:Unsigned}, rawdata, rawidx, nbytes)
    # samples are left-justified in the file
    shift = sizeof(T) - nbytes
    # first construct the native value (left-aligned)
    intval = zero(T)
    for i in 0:(nbytes-1)
        intval |= T(rawdata[rawidx+i]) << ((i+shift)*8)
    end
    intval
end

# we represent PCM values as right-justified fixed-point.
for (signed, unsigned) in [
        (Int16, UInt16),
        (Int32, UInt32),
        (Int64, UInt64)]
    @eval function decodewavbytes(::Type{Fixed{$signed, f}}, rawdata, rawidx, nbytes) where {f}
        # first decode to the left-justified unsigned int of the correct width
        intval = decodewavbytes($unsigned, rawdata, rawidx, nbytes)
        signedval = reinterpret($signed, intval)
        # f is the number of fractional bits, plus a sign bit
        nbits = f+1
        # now shift it (preserving the sign bit) so it's right-aligned
        return reinterpret(Fixed{$signed, f}, signedval >> (sizeof($unsigned)*8-nbits))
    end
end

# this doesn't do any companding, it's for PCM 8-bit-or-smaller data. WAV
# stores 8-bit PCM data in an unsigned byte ranging 0x00-0xFF, so we offset it
# to a signed fixed-point representation. If the value is smaller than 8 bits it
# will be right-aligned.
function decodewavbytes(::Type{Fixed{Int8, f}}, rawdata, rawidx, nbytes) where {f}
    @assert nbytes == 1
    nbits = f+1
    signedval = reinterpret(Int8, rawdata[rawidx]-0x80)
    return reinterpret(Fixed{Int8, f}, signedval >> (8-nbits))
end

# a bunch of other target types come from just building a unsigned integer type
# and then reinterpreting is as the target type. Use some metaprogramming
# to define all these methods. You can't reinterpret directly from a UIntX type
# into a fixed-point, so we go UInt -> Int -> Fixed
for (targettype, inttype) in (
        (Float32, UInt32),
        (Float64, UInt64),
        (Int16, UInt16),
        (Int32, UInt64),
        (Int64, UInt64))
    @eval function decodewavbytes(::Type{$targettype}, rawdata, rawidx, nbytes)
        intval = decodewavbytes($inttype, rawdata, rawidx, nbytes)
        return reinterpret($targettype, intval)
    end
end

function throw_fmt_error(fmt)
    if isextensible(fmt)
        error("Unrecognized format code (extensible): $(fmt.ext.sub_format)")
    else
        error("Unrecognized format code: $(fmt.compression_code)")
    end
end

"""
Find the format subchunk, and store up any non-data chunks it finds
along the way. This function expects the stream to be right after the end
of a previous chunk, i.e. the next think in the stream is a chunk header.
"""
function find_format(io, bytesleft)
    opt = Dict{String, Vector{Vector{UInt8}}}()
    while isnull(bytesleft) || bytesleft > SUBCHUNK_HEADER_SIZE
        (subchunk_id,
         subchunk_size,
         bytesleft) = read_subchunk_header(io, bytesleft)
        if subchunk_id == b"fmt "
            format, bytesleft = read_format(io, subchunk_size, bytesleft)
            return format, opt, bytesleft
        elseif subchunk_id == b""
            # got EOF reading the header
            break
        elseif subchunk_id == b"data"
            error("Got data chunk before fmt chunk")
        else
            chunkdata, bytesleft = read_subchunk(io, subchunk_size, bytesleft)
            pushmetadata(opt, subchunk_id, chunkdata)
        end
    end
    error("Parsed whole file without seeing a fmt chunk")
end

"""
    find_data(io, opt, bytesleft)

Find the data subchunk, and store up any non-data chunks it finds along the
way in the given `opt` dictionary. This function expects the stream to be
right after the end of a previous chunk, i.e. the next think in the stream is
a chunk header.

Returns the data subchunk size and total remaining file size, with the io
stream positioned at the beginning of the data (after the header).
"""
function find_data(io, opt, bytesleft)
    while isnull(bytesleft) || bytesleft > SUBCHUNK_HEADER_SIZE
        (subchunk_id,
         subchunk_size,
         bytesleft) = read_subchunk_header(io, bytesleft)
        if subchunk_id == b"fmt "
            warn("Got extra fmt chunk, skipping...")
            # throw away the data
            _, bytesleft = read_subchunk(io, subchunk_size, bytesleft)
        elseif subchunk_id == b"data"
            return subchunk_size, bytesleft
        elseif subchunk_id == b""
            # got EOF while reading subchunk header, bytesleft is now 0
        else
            chunkdata, bytesleft = read_subchunk(io, subchunk_size, bytesleft)
            pushmetadata(opt, subchunk_id, chunkdata)
        end
    end
    error("Parsed whole file without seeing a data chunk")
end

"""
Parse whatever is left in the WAV file and store any chunks in the given
`opt` dict.
"""
function parse_tail(io, opt, bytesleft)
    while isnull(bytesleft) || bytesleft > SUBCHUNK_HEADER_SIZE
        (subchunk_id,
         subchunk_size,
         bytesleft) = read_subchunk_header(io, bytesleft)
        if subchunk_id == b"fmt " || subchunk_id == b"data"
            warn("Got extra $(String(subchunk_id)) chunk, skipping...")
            # throw away the data
            _, bytesleft = read_subchunk(io, subchunk_size, bytesleft)
        elseif subchunk_id == b""
            # got EOF reading the header, bytesleft is now 0
        else
            chunkdata, bytesleft = read_subchunk(io, subchunk_size, bytesleft)
            pushmetadata(opt, subchunk_id, chunkdata)
        end
    end

    nothing
end

"""
    read_subchunk(io, size, bytesleft)

Read the subchunk payload with the given `size` in bytes, assuming that there
are `bytesleft` bytes left in the stream.
"""
function read_subchunk(io, size, bytesleft)
    data = read(io, size)
    if length(data) < size
        warn("Got EOF reading subchunk payload")
        bytesleft = 0
    else
        bytesleft -= size
        if isodd(size)
            # consume the padding byte
            skip(io, 1)
            bytesleft > 0 && (bytesleft -= 1)
        end
    end

    data, bytesleft
end

# The WAV specification states that numbers are written to disk in little endian form.
write_le(stream::IO, value) = write(stream, htol(value))
read_le(stream::IO, x::Type) = ltoh(read(stream, x))

"""
    read_header(io)

Read the RIFF and WAVE headers from the given `io` object
"""
function read_header(io::IO)
    # check if the given file has a valid RIFF header
    riff = Array{UInt8}(4)
    read!(io, riff)
    if riff !=  b"RIFF"
        error("Invalid RIFF file: got \"$(String(riff))\", expected \"RIFF\"")
    end

    chunk_size = Int(read_le(io, UInt32))
    if chunk_size == 0 || chunk_size == 0xffffffff
        warn("Illegal size $chunk_size in RIFF header. Reading until EOF")
        chunk_size = null
    elseif isodd(chunk_size)
        warn("Odd size $chunk_size in RIFF header, adding padding byte")
        # assume we should add a padding byte, but this is a sketchy situation
        chunk_size += 1
    end

    # check if this is a WAV file
    format = Array{UInt8}(4)
    read!(io, format)
    chunk_size -= 4
    if format != b"WAVE"
        error("Invalid WAV file: got \"$(String(format))\", expected \"WAVE\"")
    end
    # the given file size is does not include the 8-byte RIFF header, but does
    # include the WAVE string that came next
    return chunk_size
end

"""
    read_subchunk_header(io, bytesleft)

Reads a subchunk header from the given stream, assuming that there's the
given number of bytes remaining in the file. It also does some basic
validation on the subchunk.

Returns the subchunk ID, subchunk size, and new number of bytes remaining.
"""
function read_subchunk_header(io::IO, bytesleft)
    subchunk_id, subchunk_size = try
            read(io, UInt8, 4), Int(read_le(io, UInt32))
        catch ex
            ex isa EOFError || rethrow()
            warn("Got EOF expecting subchunk header")
            return b"", 0, 0
        end
    bytesleft -= 8
    if !isnull(bytesleft) && subchunk_size > bytesleft
        warn("Subchunk \"$(String(subchunk_id))\" claims to be longer than RIFF chunk length. Truncating...")
        subchunk_size = bytesleft
    end
    return subchunk_id, subchunk_size, bytesleft
end