# -*- mode: julia; -*-
Base.__precompile__(true)
module WAV
export WAVSink, WAVSource
export WAVFormatExtension, WAVFormat
export WAVE_FORMAT_PCM, WAVE_FORMAT_IEEE_FLOAT, WAVE_FORMAT_ALAW, WAVE_FORMAT_MULAW

using FileIO
import SampledSignals
using SampledSignals: SampleSource, SampleSink, nchannels
using Nulls: null, isnull

using FixedPointNumbers: Fixed

include("formats.jl")
include("companding.jl")
include("fileio.jl")

const SUBCHUNK_HEADER_SIZE = 8

"""
    WAVSink(io) # TODO: kwargs here!

A writable audio stream that writes to a WAV file represented by the given io
stream.
"""
mutable struct WAVSink{IO} <: SampleSink where IO
    io::IO
    format::WAVFormat
    totalbytes::Int
end

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
    opt::Dict{Symbol, Vector{UInt8}}
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

# SampledSignals.nframes(src::WAVSource)

function SampledSignals.unsafe_read!(src::WAVSource, buf::Array,
                                     frameoffset, framecount)
    framebytes = src.format.block_align
    samplebytes = div(framebytes, nchannels(src))
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
        src.file_bytesleft -= 0
        src.data_bytesleft -= 0
    else
        src.file_bytesleft -= bytesread
        src.data_bytesleft -= bytesread
    end
    framesread = div(bytesread, framebytes)

    copy_samples(buf, rawdata, src.format, frameoffset, framesread, samplebytes)

    if src.data_bytesleft < framebytes
        # we've reached the end of the data chunk
        if src.data_bytesleft > 0
            warn("\"data\" chunk had $(src.data_bytesleft) orphaned bytes")
            skip(src.io, src.data_bytesleft)
            src.file_bytesleft -= src.data_bytesleft
            src.data_bytesleft = 0
        end

        # grab any of the chunks that follow
        parse_tail(src.io, src.opt, src.data_bytesleft)
    end

    framesread
end

"""
    function copy_samples(dest, src,
                        format, frameoffset, framecount,
                        samplebytes)
Copy the samples from `src` (the raw UInt8 buffer from the stream) into `dest`
(the array used for unsafe_read)
"""
function copy_samples(dest::Array{T, N}, src,
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
function copy_samples(dest::Array{Fixed{Int16, 15}, N}, src,
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

"""
Convert the underlying wav format to a Julia type that we can
use to parameterize sources/sinks for type-stable reads and writes.

Companded streams are represented as 16-bit fixed-point because that's the
user-facing representation.
"""
function type_from_format(fmt)
    nbits = bits_per_sample(fmt)
    if isformat(fmt, WAVE_FORMAT_PCM)
        if nbits <= 8
            return Fixed{Int8, nbits-1}
        elseif nbits <= 16
            return Fixed{Int16, nbits-1}
        elseif nbits <= 32
            return Fixed{Int32, nbits-1}
        elseif nbits <= 64
            return Fixed{Int64, nbits-1}
        else
            error("Only support PCM up to 64-bit")
        end
    elseif isformat(fmt, WAVE_FORMAT_IEEE_FLOAT)
        if nbits == 32
            return Float32
        elseif nbits == 64
            return Float64
        else
            error("FORMAT_IEEE_FLOAT only supports 32 and 64-bit floats")
        end
    # Companded streams will appear as a 16-bit fixed-point stream, and we'll
    # compand to/from 8-bit on read/write. The spec says that companded streams
    # are always 8-bits/sample
    elseif isformat(fmt, WAVE_FORMAT_MULAW)
        return Fixed{Int16, 15}
    elseif isformat(fmt, WAVE_FORMAT_ALAW)
        return Fixed{Int16, 15}
    else
        throw_fmt_error(fmt)
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
    opt = Dict{Symbol, Vector{UInt8}}()
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
            opt[Symbol(subchunk_id)], bytesleft = read_subchunk(io,
                                                                subchunk_size,
                                                                bytesleft)
        end
    end
    error("Parsed whole file without seeing a fmt chunk")
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
            opt[Symbol(subchunk_id)], bytesleft = read_subchunk(io,
                                                                subchunk_size,
                                                                bytesleft)
        end
    end

    nothing
end

"""
Find the next data subchunk, and store up any non-data chunks it finds along the
way in the given `opt` dictionary. This function expects the stream to be right
after the end of a previous chunk, i.e. the next think in the stream is a chunk
header.

Returns the data subchunk size and total remaining file size, with the io stream
positioned at the beginning of the data (after the header).
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
            opt[Symbol(subchunk_id)], bytesleft = read_subchunk(io,
                                                                subchunk_size,
                                                                bytesleft)
        end
    end
    error("Parsed whole file without seeing a data chunk")
end

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
Reads a chunk header from the given stream, assuming that there's the
given number of bytes remaining in the file. It also does some basic validation
on the subchunk.

Returns the subchunk ID, subchunk size, and new number of bytes remaining.
"""
function read_subchunk_header(io::IO, bytesleft)
    subchunk_id, subchunk_size = try
            read(io, UInt8, 4), Int(read_le(io, UInt32))
        catch ex
            ex isa EOFError || rethrow()
            warn("Got EOF reading subchunk header")
            return b"", 0, 0
        end
    bytesleft -= 8
    if !isnull(bytesleft) && subchunk_size > bytesleft
        warn("Subchunk \"$(String(subchunk_id))\" claims to be longer than RIFF chunk length. Truncating...")
        subchunk_size = bytesleft
    end
    return subchunk_id, subchunk_size, bytesleft
end

function write_header(io::IO, databytes)
    write(io, b"RIFF") # RIFF header
    write_le(io, databytes) # chunk_size
    write(io, b"WAVE")
end
write_standard_header(io, databytes) = write_header(io, UInt32(databytes + 36))
write_extended_header(io, databytes) = write_header(io, UInt32(databytes + 60))

function ieee_float_container_type(nbits)
    if nbits == 32
        return Float32
    elseif nbits == 64
        return Float64
    end

    error("$nbits bits is not supported for WAVE_FORMAT_IEEE_FLOAT.")
end

function write_pcm_samples(io::IO, fmt::WAVFormat, samples::AbstractArray{<:Integer})
    const nbits = bits_per_sample(fmt)
    # number of bytes per sample
    nbytes = ceil(Integer, nbits / 8)
    for i = 1:size(samples, 1)
        for j = 1:size(samples, 2)
            my_sample = samples[i, j]
            # shift my_sample into the N most significant bits
            my_sample <<= nbytes * 8 - nbits
            for k = 1:nbytes
                write_le(io, convert(UInt8, my_sample & 0xff))
                my_sample = my_sample >> 8
            end
        end
    end
end

function write_pcm_samples(io::IO, fmt::WAVFormat, samples::AbstractArray{T}) where T <: AbstractFloat
    nbits = bits_per_sample(fmt)
    # Scale the floating point values to the PCM range
    if nbits > 8
        # two's complement
        samples = convert(Array{pcm_container_type(nbits)}, round.(samples * (2.0^(nbits - 1) - 1)))
    else
        # offset binary
        samples = convert(Array{UInt8}, round.((samples .+ 1.0) / 2.0 * (2.0^nbits - 1)))
    end
    return write_pcm_samples(io, fmt, samples)
end

function write_ieee_float_samples(io::IO, samples)
    # Interleave the channel samples before writing to the stream.
    for i = 1:size(samples, 1) # for each sample
        for j = 1:size(samples, 2) # for each channel
            write_le(io, samples[i, j])
        end
    end
end

# take the loop variable type out of the loop
function write_ieee_float_samples(io::IO, fmt::WAVFormat, samples)
    floatType = ieee_float_container_type(bits_per_sample(fmt))
    write_ieee_float_samples(io, convert(Array{floatType}, samples))
end

function write_data(io::IO, fmt::WAVFormat, samples::AbstractArray)
    if isformat(fmt, WAVE_FORMAT_PCM)
        return write_pcm_samples(io, fmt, samples)
    elseif isformat(fmt, WAVE_FORMAT_IEEE_FLOAT)
        return write_ieee_float_samples(io, fmt, samples)
    elseif isformat(fmt, WAVE_FORMAT_MULAW)
        return write_companded_samples(io, samples, compress_sample_mulaw)
    elseif isformat(fmt, WAVE_FORMAT_ALAW)
        return write_companded_samples(io, samples, compress_sample_alaw)
    else
        error("$(fmt.compression_code) is an unsupported compression code.")
    end
end

make_range(subrange) = subrange
make_range(subrange::Number) = 1:convert(Int, subrange)

function wavread(filename::AbstractString; subrange=Void, format="double")
    open(filename, "r") do io
        wavread(io, subrange=subrange, format=format)
    end
end

# These are the MATLAB compatible signatures
wavread(filename::AbstractString, fmt::AbstractString) = wavread(filename, format=fmt)
wavread(filename::AbstractString, n) = wavread(filename, subrange=n)
wavread(filename::AbstractString, n, fmt) = wavread(filename, subrange=n, format=fmt)

get_default_compression(::AbstractArray{T}) where T <: Integer = WAVE_FORMAT_PCM
get_default_compression(::AbstractArray{T}) where T<:AbstractFloat = WAVE_FORMAT_IEEE_FLOAT
get_default_pcm_precision(::AbstractArray{UInt8}) = 8
get_default_pcm_precision(::AbstractArray{Int16}) = 16
get_default_pcm_precision(::Any) = 24

function get_default_precision(samples, compression)
    if compression == WAVE_FORMAT_ALAW || compression == WAVE_FORMAT_MULAW
        return 8
    elseif compression == WAVE_FORMAT_IEEE_FLOAT
        return 32
    end
    get_default_pcm_precision(samples)
end

function wavwrite(samples::AbstractArray, io::IO; Fs=8000, nbits=0, compression=0,
                  chunks::Dict{Symbol, Array{UInt8,1}}=Dict{Symbol, Array{UInt8,1}}())
    if compression == 0
        compression = get_default_compression(samples)
    elseif compression == WAVE_FORMAT_ALAW || compression == WAVE_FORMAT_MULAW
        nbits = 8
    end
    if nbits == 0
        nbits = get_default_precision(samples, compression)
    end
    compression_code = compression
    nchannels = size(samples, 2)
    sample_rate = Fs
    my_nbits = ceil(Integer, nbits / 8) * 8
    block_align = my_nbits / 8 * nchannels
    bps = sample_rate * block_align
    data_length::UInt32 = size(samples, 1) * block_align
    ext = WAVFormatExtension()

    if nchannels > 2 || my_nbits > 16 || my_nbits != nbits
        compression_code = WAVE_FORMAT_EXTENSIBLE
        valid_bits_per_sample = nbits
        channel_mask = 0
        sub_format = Array{UInt8, 1}(0)
        if compression == WAVE_FORMAT_PCM
            sub_format = KSDATAFORMAT_SUBTYPE_PCM
        elseif compression == WAVE_FORMAT_IEEE_FLOAT
            sub_format = KSDATAFORMAT_SUBTYPE_IEEE_FLOAT
        elseif compression == WAVE_FORMAT_ALAW
            sub_format = KSDATAFORMAT_SUBTYPE_ALAW
        elseif compression == WAVE_FORMAT_MULAW
            sub_format = KSDATAFORMAT_SUBTYPE_MULAW
        else
            error("Unsupported extension sub format: $compression")
        end
        ext = WAVFormatExtension(valid_bits_per_sample, channel_mask, sub_format)
        write_extended_header(io, data_length)
    else
        write_standard_header(io, data_length)
    end
    fmt = WAVFormat(compression_code,
                    nchannels,
                    sample_rate,
                    bps,
                    block_align,
                    my_nbits,
                    ext)
    write_format(io, fmt)

    for eachchunk in chunks
        write(io, eachchunk[1])
        write_le(io, UInt32(length(eachchunk[2])))
        for eachbyte in eachchunk[2]
            write(io, eachbyte)
        end
    end

    # write the data subchunk header
    write(io, b"data")
    write_le(io, data_length) # UInt32
    write_data(io, fmt, samples)
end

function wavwrite(samples::AbstractArray, filename::AbstractString; Fs=8000, nbits=0, compression=0,
                  chunks::Dict{Symbol, Array{UInt8,1}}=Dict{Symbol, Array{UInt8,1}}())
    open(filename, "w") do io
        wavwrite(samples, io, Fs=Fs, nbits=nbits, compression=compression, chunks=chunks)
    end
end

# function wavappend(samples::AbstractArray, io::IO)
#     seekstart(io)
#     chunk_size = read_header(io)
#     subchunk_id = Array{UInt8}(4)
#     read!(io, subchunk_id)
#     subchunk_size = read_le(io, UInt32)
#     if subchunk_id != b"fmt "
#         error("First chunk is not the format")
#     end
#     fmt = read_format(io, subchunk_size)

#     if fmt.nchannels != size(samples,2)
#         error("Number of channels do not match")
#     end

#     data_length = size(samples, 1) * fmt.block_align

#     seek(io,4)
#     write_le(io, convert(UInt32, chunk_size + data_length))

#     seek(io,64)
#     subchunk_size = read_le(io, UInt32)
#     seek(io,64)
#     write_le(io, convert(UInt32, subchunk_size + data_length))

#     seekend(io)
#     write_data(io, fmt, samples)
# end

# function wavappend(samples::AbstractArray, filename::AbstractString)
#     open(filename, "a+") do io
#         wavappend(samples,io)
#     end
# end

wavwrite(y::AbstractArray, f::Real, filename::AbstractString) = wavwrite(y, filename, Fs=f)
wavwrite(y::AbstractArray, f::Real, n::Real, filename::AbstractString) = wavwrite(y, filename, Fs=f, nbits=n)

# support for writing native arrays...
wavwrite(y::AbstractArray{T}, io::IO) where T <: Integer = wavwrite(y, io, nbits=sizeof(T)*8)
wavwrite(y::AbstractArray{T}, filename::AbstractString) where T <: Integer = wavwrite(y, filename, nbits=sizeof(T)*8)
wavwrite(y::AbstractArray{Int32}, io::IO) = wavwrite(y, io, nbits=24)
wavwrite(y::AbstractArray{Int32}, filename::AbstractString) = wavwrite(y, filename, nbits=24)
wavwrite(y::AbstractArray{T}, io::IO) where T <: AbstractFloat = wavwrite(y, io, nbits=sizeof(T)*8, compression=WAVE_FORMAT_IEEE_FLOAT)
wavwrite(y::AbstractArray{T}, filename::AbstractString) where T <: AbstractFloat = wavwrite(y, filename, nbits=sizeof(T)*8, compression=WAVE_FORMAT_IEEE_FLOAT)

end # module
