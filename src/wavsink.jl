"""
    function WAVSink(io; eltype=nothing,
                         nchannels=2,
                         samplerate=48000,
                         compression=:none)

A writable audio stream that writes to a WAV file represented by the given io
stream.

`compression` should be one of `:none`, `:mulaw`, or `:alaw`, and if not
`:none` will write a wav file with 8-bit samples encoded with the given
compression scheme.

`eltype` will default to `Float32` unless `compression` is used, in which
case it will be `Fixed{Int16, 15}`, i.e. it will expect 16-bit signed
fixed-point samples. For floating-point and PCM samples you can specify a
different floating or fixed-point sample type, respectively.

## Examples:

```julia
t = linspace(0, 2, 16001)
sig = cos.(2pi*220t)
open("/tmp/outfile.wav") do io
    WAVSink(io; samplerate=8000, nchannels=1, compression=:mulaw) do sink
        write(sink, sig)
    end
end
```

You can also provide a filename rather than an `IO` object, and `WAVSink` will
handle opening and closing the file:

```julia
t = linspace(0, 2, 16001)
sig = cos.(2pi*220t)
WAVSink("/tmp/outfile.wav") do sink
    write(sink, sig)
end
```

Most people will interact with `WAVSink` using FileIO's `save` function instead,
which takes the same keyword arguments. When using `save` the properties of the
given signal will override the `WAVSink` defaults listed above, so the data will
be saved with minimal information loss. If the signal is given as an `Array`
then the element type and channel count will come from the array, with a default
sampling rate. If a `SampleBuf` (from SampledSignals.jl) is used, the file will
also match the buffer's sampling rate.

```julia
t = linspace(0, 2, 16001)
sig = cos.(2pi*220t)
save("/tmp/outfile.wav", sig; samplerate=8000)
```
"""
mutable struct WAVSink{IO} <: SampleSink where IO
    io::IO
    format::WAVFormat
    riffpos::UInt64 # position of the RIFF header length field
    rifflen::UInt32 # total bytes written within RIFF chunk
    datapos::UInt64 # position of the data header length field
    datalen::UInt32 # total bytes written within data chunk
    factpos::UInt64 # position of the fact header length field
    factframes::UInt32 # total frames written within data chunk (for fact chunk)
    writingdata::Bool # true if we're in the middle of writing the data chunk
    # store up metadata chunks that will get written before we close the file
    pending_metadata::Vector{Tuple{ChunkKey, Vector{UInt8}}}
end

# TODO: support user-specified channel mapping
# TODO: allow length hints so we can provide correct lengths at the beginning
# if we know them

# use DUMMY_LENGTH when we need to write a header but don't know how long the
# content will be
const DUMMY_LENGTH = 0xffffffff

function WAVSink(io::T; eltype=nothing,
                        nchannels=2,
                        samplerate=48000,
                        compression=:none) where {T <: IO}
    if !(compression in COMPRESSIONS)
        error("`compression` must be one of $COMPRESSIONS")
    end
    if eltype == nothing
        if compression != :none
            eltype = Fixed{Int16, 15}
        else
            eltype = Float32
        end
    end
    if compression != :none && eltype != Fixed{Int16, 15}
        error("a-law and mu-law compression is only compatible with Fixed{Int16, 16} samples. Leave `eltype` unspecified for default.")
    end
    riffpos = write_header(io, DUMMY_LENGTH)
    rifflen = 4 # needs to include the initial "WAVE" block
    fmt = WAVFormat(eltype, nchannels, samplerate, compression)
    rifflen += write_format(io, fmt)

    if fmt != WAVE_FORMAT_PCM
        # non-PCM files require a FACT chunk
        write(io, b"fact")
        write_le(io, UInt32(4))
        factpos = position(io)
        write_le(io, UInt32(0xffffffff))
        rifflen += 12
    else
        factpos = 0
    end
    datapos = 0
    datalen = 0
    factframes = 0
    writingdata = false
    pendingchunks = Tuple{ChunkKey, Vector{UInt8}}[]

    # now we're ready to write more chunks from the user
    WAVSink{T}(io, fmt, riffpos, rifflen,
            datapos, datalen, factpos, factframes,
            writingdata, Tuple{ChunkKey, Vector{UInt8}}[])
end

SampledSignals.nchannels(sink::WAVSink) = Int(sink.format.nchannels)
SampledSignals.eltype(sink::WAVSink) = type_from_format(sink.format)
SampledSignals.samplerate(sink::WAVSink) = Int(sink.format.sample_rate)

function write_dataheader(sink)
    write(sink.io, b"data")
    sink.datapos = position(sink.io)
    write_le(sink.io, UInt32(0xffffffff))
    sink.rifflen += 8
end

# implement SampledSignals's write method to plug into the stream conversion
# infrastructure
function SampledSignals.unsafe_write(sink::WAVSink, buf::Array,
                                     frameoffset, framecount)
    if !sink.writingdata
        write_dataheader(sink)
        sink.writingdata = true
    end

    framebytes = sink.format.block_align
    samplebytes = framebytes ÷ nchannels(sink)
    bytestowrite = framecount * framebytes
    # TODO: pre-allocate this buffer
    rawdata = Vector{UInt8}(bytestowrite)

    encode_buffer(rawdata, buf, sink.format, frameoffset, framecount, samplebytes)
    write(sink.io, rawdata)
    sink.rifflen += bytestowrite
    sink.datalen += bytestowrite
    sink.factframes += framecount
end

"""
    writemeta(sink::WAVSink, key::ChunkKey, chunkdata::Vector{UInt8})

Write a metadata chunk to the given WAV file. If no sample data has been
written yet than this chunk will be written immediately (and will precede the
sample data in the resulting file). If sample data was previously written to
the stream, this chunk will be queued and written after the sample data when
the stream is closed.
"""
function writemeta(sink::WAVSink, key::ChunkKey, chunkdata::Vector{UInt8})
    if sink.writingdata
        # we're in the middle of the sample data, so defer this metadata
        # to the end
        push!(sink, (key, chunkdata))
    else
        sink.rifflen += writechunk(sink.io, key, chunkdata)
    end
end

hasfact(sink::WAVSink) = sink.factpos != 0
hasdata(sink::WAVSink) = sink.datapos != 0

function Base.close(sink::WAVSink)
    if !sink.writingdata
        # looks like we're closing without having written any data, so we'll
        # need to write the data header here. The length will get set later
        write_dataheader(sink)
    end
    if isodd(sink.datalen)
        # add the padding byte for the sample data
        write(sink.io, 0x00)
        sink.rifflen += 1
    end
    for (key, chunkdata) in sink.pending_metadata
        sink.rifflen += writechunk(sink.io, key, chunkdata)
    end
    # there's currently no way to check to see if the stream is actually
    # seekable, so we just try and abort if it doesn't work out.
    try
        currentpos = position(sink.io)
        seek(sink.io, sink.riffpos)
        write_le(sink.io, sink.rifflen)
        if hasfact(sink)
            seek(sink.io, sink.factpos)
            write_le(sink.io, sink.factframes)
        end
        if hasdata(sink)
            seek(sink.io, sink.datapos)
            write_le(sink.io, sink.datalen)
        end
        seek(sink.io, currentpos)
    catch ex
        warn("Caught error trying to seek back to write length fields: $ex")
    end

    nothing
end

"""
    function encode_buffer(dest, src,
                           format, frameoffset, framecount,
                           samplebytes)
Copy the samples from `src` (the buffer of samples from `unsafe_write`) into
`dest` (the buffer of `UInt8` that we'll write to the file)
"""
function encode_buffer(dest, src::Array{T, N},
                       format, frameoffset, framecount,
                       samplebytes) where {T, N}
    rawidx = 1
    for frame in (1:framecount) + frameoffset
        for ch in 1:size(src, 2)
            encodewavbytes(src[frame, ch], dest, rawidx, samplebytes)
            rawidx += samplebytes
        end
    end
end

# we need special handling for writing from 16-bit because we need to
# check whether we actually have an 8-bit companded stream
function encode_buffer(dest, src::Array{Fixed{Int16, 15}, N},
                       format, frameoffset, framecount,
                       samplebytes) where N
    rawidx = 1
    if isformat(format, WAVE_FORMAT_PCM)
        # OK, just a normal 16-bit PCM
        for frame in (1:framecount) + frameoffset, ch in 1:size(src, 2)
            encodewavbytes(src[frame, ch], dest, rawidx, 2)
            rawidx += 2
        end
    elseif isformat(format, WAVE_FORMAT_MULAW)
        for frame in (1:framecount) + frameoffset, ch in 1:size(src, 2)
            dest[rawidx] = encodemulaw(src[frame, ch])
            rawidx += 1
        end
    elseif isformat(format, WAVE_FORMAT_ALAW)
        for frame in (1:framecount) + frameoffset, ch in 1:size(src, 2)
            dest[rawidx] = encodealaw(src[frame, ch])
            rawidx += 1
        end
    else
        throw_fmt_error(format)
    end

    framecount
end
# TODO: make sure we test odd-sized chunks
# internal function that writes the given chunk directly to the IO stream
function writechunk(io::IO, key::Vector{UInt8}, chunkdata::Vector{UInt8})
    nwritten = 0
    @assert length(key) == 4
    write(io, key)
    write_le(io, UInt32(length(chunkdata)))
    nwritten += 8
    write(io, chunkdata)
    nwritten += length(chunkdata)
    if isodd(length(chunkdata))
        # add the padding byte
        write(io, UInt8(0))
        nwritten += 1
    end

    nwritten
end

# write the header block and return the stream position of the RIFF length
# field so we can seek back to it later if we need
function write_header(io::IO, databytes)
    write(io, b"RIFF") # RIFF header
    lenpos = position(io)
    write_le(io, databytes) # chunk_size
    write(io, b"WAVE")

    lenpos
end

"""
    encodewavbytes(value, rawdata, rawidx, nbytes)

Encode the given `value` into raw bytes to be written to a wav file. The data
is written to the `rawdata` vector, starting at `rawidx`. the `nbytes`
argument specifies how many bytes the sample takes up in the buffer, (e.g.
24-bit audio takes 3 bytes).

These functions should be type-stable and fast, so they can be called in an
inner loop.
"""
function encodewavbytes end

# write down the least-significant `nbytes` bytes in little-endian order
function encodewavbytes(value::T, rawdata, rawidx, nbytes) where {T<:Integer}
    for i in 0:(nbytes-1)
        rawdata[rawidx+i] = (value >> 8i) & 0xff
    end
end

# these types are all encoded by just reinterpreting to a different type and
# then encoding that type
for (inttype, fromtype) in [
        (UInt64, Float64),
        (UInt32, Float32)]
    @eval function encodewavbytes(value::$fromtype, rawdata, rawidx, nbytes)
        encodewavbytes(reinterpret($inttype, value), rawdata, rawidx, nbytes)
    end
end

function encodewavbytes(value::Fixed{T, f}, rawdata, rawidx, nbytes) where {T, f}
    intval = reinterpret(T, value)
    # fixed-point values should be left-justified within the byte
    shift = 8nbytes - (f+1)
    encodewavbytes(intval << shift, rawdata, rawidx, nbytes)
end

# 8-bit(and less) PCM is treated as unsigned 0x00 - 0xff, so we offset it so to
# users it looks like a normal 8-bit fixed-point value
function encodewavbytes(value::Fixed{Int8, f}, rawdata, rawidx, nbytes) where {f}
    @assert nbytes == 1
    intval = reinterpret(UInt8, reinterpret(Int8, value))
    # fixed-point values should be left-justified within the byte
    shift = 8 - (f+1)
    rawdata[rawidx] = intval << shift + 0x80
end

# function wavread(filename::AbstractString; subrange=Void, format="double")
#     open(filename, "r") do io
#         wavread(io, subrange=subrange, format=format)
#     end
# end
#
# get_default_compression(::AbstractArray{T}) where T <: Integer = WAVE_FORMAT_PCM
# get_default_compression(::AbstractArray{T}) where T<:AbstractFloat = WAVE_FORMAT_IEEE_FLOAT
# get_default_pcm_precision(::AbstractArray{UInt8}) = 8
# get_default_pcm_precision(::AbstractArray{Int16}) = 16
# get_default_pcm_precision(::Any) = 24
#
# function get_default_precision(samples, compression)
#     if compression == WAVE_FORMAT_ALAW || compression == WAVE_FORMAT_MULAW
#         return 8
#     elseif compression == WAVE_FORMAT_IEEE_FLOAT
#         return 32
#     end
#     get_default_pcm_precision(samples)
# end
#
# function wavwrite(samples::AbstractArray, io::IO; Fs=8000, nbits=0, compression=0,
#                   chunks::Dict{Symbol, Array{UInt8,1}}=Dict{Symbol, Array{UInt8,1}}())
#     if compression == 0
#         compression = get_default_compression(samples)
#     elseif compression == WAVE_FORMAT_ALAW || compression == WAVE_FORMAT_MULAW
#         nbits = 8
#     end
#     if nbits == 0
#         nbits = get_default_precision(samples, compression)
#     end
#     compression_code = compression
#     nchannels = size(samples, 2)
#     sample_rate = Fs
#     my_nbits = ceil(Integer, nbits / 8) * 8
#     block_align = my_nbits / 8 * nchannels
#     bps = sample_rate * block_align
#     data_length::UInt32 = size(samples, 1) * block_align
#     ext = WAVFormatExtension()
#
#     if nchannels > 2 || my_nbits > 16 || my_nbits != nbits
#         compression_code = WAVE_FORMAT_EXTENSIBLE
#         valid_bits_per_sample = nbits
#         channel_mask = 0
#         sub_format = Array{UInt8, 1}(0)
#         if compression == WAVE_FORMAT_PCM
#             sub_format = KSDATAFORMAT_SUBTYPE_PCM
#         elseif compression == WAVE_FORMAT_IEEE_FLOAT
#             sub_format = KSDATAFORMAT_SUBTYPE_IEEE_FLOAT
#         elseif compression == WAVE_FORMAT_ALAW
#             sub_format = KSDATAFORMAT_SUBTYPE_ALAW
#         elseif compression == WAVE_FORMAT_MULAW
#             sub_format = KSDATAFORMAT_SUBTYPE_MULAW
#         else
#             error("Unsupported extension sub format: $compression")
#         end
#         ext = WAVFormatExtension(valid_bits_per_sample, channel_mask, sub_format)
#         write_extended_header(io, data_length)
#     else
#         write_standard_header(io, data_length)
#     end
#     fmt = WAVFormat(compression_code,
#                     nchannels,
#                     sample_rate,
#                     bps,
#                     block_align,
#                     my_nbits,
#                     ext)
#     write_format(io, fmt)
#
#     for eachchunk in chunks
#         write(io, eachchunk[1])
#         write_le(io, UInt32(length(eachchunk[2])))
#         for eachbyte in eachchunk[2]
#             write(io, eachbyte)
#         end
#     end
#
#     # write the data subchunk header
#     write(io, b"data")
#     write_le(io, data_length) # UInt32
#     write_data(io, fmt, samples)
# end
#
# function wavwrite(samples::AbstractArray, filename::AbstractString; Fs=8000, nbits=0, compression=0,
#                   chunks::Dict{Symbol, Array{UInt8,1}}=Dict{Symbol, Array{UInt8,1}}())
#     open(filename, "w") do io
#         wavwrite(samples, io, Fs=Fs, nbits=nbits, compression=compression, chunks=chunks)
#     end
# end

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
