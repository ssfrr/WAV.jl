# -*- mode: julia; -*-
Base.__precompile__(true)
module WAV
export WAVFormatExtension, WAVFormat
export WAVE_FORMAT_PCM, WAVE_FORMAT_IEEE_FLOAT, WAVE_FORMAT_ALAW, WAVE_FORMAT_MULAW
using FileIO

include("companding.jl")
include("formats.jl")

mutable struct WAVSink{IO} <: SampleSink where IO
    io::IO
end

mutable struct WAVSource{IO} <: SampleSource where IO
    io::IO
    format::WAVFormat
    opt::Dict{Symbol, Vector{UInt8}}
    bytesleft::Int
end

function WAVSource(io)
    bytesleft = read_wave_header(io)
    fmt, opts, bytesleft = find_format(io, bytesleft)
    eltype = type_from_format(fmt)

    # now we have a WAVSource ready to provide audio data
    src
end

"""
Convert the underlying wav format to a Julia type that we can
use to parameterize sources/sinks for type-stable reads and writes.
"""
function type_from_format(fmt)
    if isformat(fmt, WAVE_FORMAT_PCM)
        return Fixed{Int16, 15}
    elseif isformat(fmt, WAVE_FORMAT_IEEE_FLOAT)
        return Float32
    elseif isformat(fmt, WAVE_FORMAT_MULAW)
        # This stream will appear as a 16-bit fixed-point stream,
        # and we'll compand to/from 8-bit on read/write
        return Fixed{Int16, 15}
    elseif isformat(fmt, WAVE_FORMAT_ALAW)
        return Fixed{Int16, 15}
    else
        if isextensible(fmt)
            error("Unrecognized format code (extensible): $(fmt.ext.sub_format)")
        else
            error("Unrecognized format code: $(fmt.compression_code)")
        end
    end
end

"""
Find the format subchunk, and store up any non-data chunks it finds
along the way.
"""
function find_format(io, bytesleft)
    opt = Dict{Symbol, Vector{UInt8}}()
    while bytesleft > 0
        (subchunk_id,
         subchunk_size,
         bytesleft) = read_subchunk_header(io, bytesleft)
        if subchunk_id == b"fmt "
            format = read_format(io, subchunk_size)
            bytesleft -= subchunk_size
            return format, opt, bytesleft
        elseif subchunk_id == b"data"
            error("Got data chunk before fmt chunk")
        elseif subchunk_id == b""
            # corrupt chunk, don't do anything
        else
            opt[Symbol(subchunk_id)] = read(io, UInt8, subchunk_size)
            bytesleft -= subchunk_size
        end
    end
    error("Parsed whole file without seeing a fmt chunk")
end

wavplay(fname) = wavplay(wavread(fname)[1:2]...)

# The WAV specification states that numbers are written to disk in little endian form.
write_le(stream::IO, value) = write(stream, htol(value))
read_le(stream::IO, x::Type) = ltoh(read(stream, x))

function read_wave_header(io::IO)
    # check if the given file has a valid RIFF header
    riff = Array{UInt8}(4)
    read!(io, riff)
    if riff !=  b"RIFF"
        error("Invalid WAV file: The RIFF header is invalid")
    end

    chunk_size = Int(read_le(io, UInt32))

    # check if this is a WAV file
    format = Array{UInt8}(4)
    read!(io, format)
    if format != b"WAVE"
        error("Invalid WAV file: the format is not WAVE")
    end
    # the given file size doesn't include the "RIFF" and "WAVE" strings,
    # but does include the 32-bit file size field itself, so we want
    # to compensate for that
    return chunk_size - 4
end

"""
Reads a chunk header from the given stream, assuming that there's the
given number of bytes remaining in the file.

Returns the subchunk ID, subchunk size, and new number of bytes remaining.
"""
function read_subchunk_header(io::IO, bytesleft)
    if bytesleft < 8
        warn("File ended with partial subchunk header")
        return b"", 0, 0
    end
    subchunk_id = read(io, UInt8, 4)
    subchunk_size = Int(read_le(io, UInt32))
    bytesleft -= 8
    if subchunk_size > bytesleft
        warn("File ended with partial subchunk data")
        subchunk_size = bytesleft
    end
    return subchunk_id, subchunk_size, bytesleft
end

function write_header(io::IO, data_length::UInt32)
    write(io, b"RIFF") # RIFF header
    write_le(io, data_length) # chunk_size
    write(io, b"WAVE")
end
write_standard_header(io, data_length) = write_header(io, UInt32(data_length + 36))
write_extended_header(io, data_length) = write_header(io, UInt32(data_length + 60))

function read_format(io::IO, chunk_size::UInt32)
    # can I read in all of the fields at once?
    orig_chunk_size = convert(Int, chunk_size)
    if chunk_size < 16
        error("The WAVE Format chunk must be at least 16 bytes")
    end
    compression_code = read_le(io, UInt16)
    nchannels = read_le(io, UInt16)
    sample_rate = read_le(io, UInt32)
    bytes_per_second = read_le(io, UInt32)
    block_align = read_le(io, UInt16)
    nbits = read_le(io, UInt16)
    ext = Array{UInt8, 1}(0)
    chunk_size -= 16
    if chunk_size > 0
        extra_bytes_length = read_le(io, UInt16)
        if extra_bytes_length == 22
            ext = Array{UInt8}(extra_bytes_length)
            read!(io, ext)
        end
    end
    return WAVFormat(compression_code,
                     nchannels,
                     sample_rate,
                     bytes_per_second,
                     block_align,
                     nbits,
                     WAVFormatExtension(ext))
end

function write_format(io::IO, fmt::WAVFormat)
    len = 16 # 16 is size of base format chunk
    if isextensible(fmt)
        len += 24 # 24 is the added length needed to encode the extension
    end
    # write the fmt subchunk header
    write(io, b"fmt ")
    write_le(io, convert(UInt32, len)) # subchunk length

    write_le(io, fmt.compression_code) # audio format (UInt16)
    write_le(io, fmt.nchannels) # number of channels (UInt16)
    write_le(io, fmt.sample_rate) # sample rate (UInt32)
    write_le(io, fmt.bytes_per_second) # byte rate (UInt32)
    write_le(io, fmt.block_align) # byte align (UInt16)
    write_le(io, fmt.nbits) # number of bits per sample (UInt16)

    if isextensible(fmt)
        write_le(io, convert(UInt16, 22))
        write_le(io, fmt.ext.nbits)
        write_le(io, fmt.ext.channel_mask)
        @assert length(fmt.ext.sub_format) == 16
        write(io, fmt.ext.sub_format)
    end
end

function pcm_container_type(nbits::Unsigned)
    if nbits > 32
        return Int64
    elseif nbits > 16
        return Int32
    elseif nbits > 8
        return Int16
    end
    return  UInt8
end

ieee_float_container_type(nbits) = (nbits == 32 ? Float32 : (nbits == 64 ? Float64 : error("$nbits bits is not supported for WAVE_FORMAT_IEEE_FLOAT.")))

function read_pcm_samples(io::IO, fmt::WAVFormat, subrange)
    nbits = bits_per_sample(fmt)
    if isempty(subrange)
        return Array{pcm_container_type(nbits), 2}(0, fmt.nchannels)
    end
    samples = Array{pcm_container_type(nbits), 2}(length(subrange), fmt.nchannels)
    sample_type = eltype(samples)
    nbytes = ceil(Integer, nbits / 8)
    bitshift = [0x0, 0x8, 0x10, 0x18, 0x20, 0x28, 0x30, 0x38, 0x40]
    mask = UInt64(0x1) << (nbits - 1)
    if nbits <= 8
        mask = UInt64(0)
    end
    skip(io, convert(UInt, (first(subrange) - 1) * nbytes * fmt.nchannels))
    for i = 1:size(samples, 1)
        for j = 1:size(samples, 2)
            raw_sample = Array{UInt8}(nbytes)
            read!(io, raw_sample)
            my_sample = UInt64(0)
            for k = 1:nbytes
                my_sample |= convert(UInt64, raw_sample[k]) << bitshift[k]
            end
            my_sample >>= nbytes * 8 - nbits
            # sign extend negative values
            my_sample = xor(my_sample, mask) - mask
            samples[i, j] = convert(sample_type, signed(my_sample))
        end
    end
    samples
end

function read_ieee_float_samples(io::IO, fmt::WAVFormat, subrange, floatType)
    if isempty(subrange)
        return Array{floatType, 2}(0, fmt.nchannels)
    end
    nblocks = length(subrange)
    samples = Array{floatType, 2}(nblocks, fmt.nchannels)
    nbits = bits_per_sample(fmt)
    skip(io, convert(UInt, (first(subrange) - 1) * (nbits / 8) * fmt.nchannels))
    for i = 1:nblocks
        for j = 1:fmt.nchannels
            samples[i, j] = read_le(io, floatType)
        end
    end
    samples
end

# take the loop variable type out of the loop
function read_ieee_float_samples(io::IO, fmt::WAVFormat, subrange)
    floatType = ieee_float_container_type(bits_per_sample(fmt))
    read_ieee_float_samples(io, fmt, subrange, floatType)
end

# PCM data is two's-complement except for resolutions of 1-8 bits, which are represented as offset binary.

# support every bit width from 1 to 8 bits
convert_pcm_to_double(samples::AbstractArray{UInt8}, nbits::Integer) = convert(Array{Float64}, samples) ./ (2.0^nbits - 1) .* 2.0 .- 1.0
convert_pcm_to_double(::AbstractArray{Int8}, ::Integer) = error("WAV files use offset binary for less than 9 bits")
# support every bit width from 9 to 64 bits
convert_pcm_to_double(samples::AbstractArray{T}, nbits::Integer) where T <: Signed = convert(Array{Float64}, samples) / (2.0^(nbits - 1) - 1)

function read_data(io::IO, chunk_size, fmt::WAVFormat, format, subrange)
    # "format" is the format of values, while "fmt" is the WAV file level format
    convert_to_double = x -> convert(Array{Float64}, x)

    if subrange === Void
        # each block stores fmt.nchannels channels
        subrange = 1:convert(UInt, chunk_size / fmt.block_align)
    end
    if isformat(fmt, WAVE_FORMAT_PCM)
        samples = read_pcm_samples(io, fmt, subrange)
        convert_to_double = x -> convert_pcm_to_double(x, bits_per_sample(fmt))
    elseif isformat(fmt, WAVE_FORMAT_IEEE_FLOAT)
        samples = read_ieee_float_samples(io, fmt, subrange)
    elseif isformat(fmt, WAVE_FORMAT_MULAW)
        samples = read_mulaw_samples(io, fmt, subrange)
        convert_to_double = x -> convert_pcm_to_double(x, 16)
    elseif isformat(fmt, WAVE_FORMAT_ALAW)
        samples = read_alaw_samples(io, fmt, subrange)
        convert_to_double = x -> convert_pcm_to_double(x, 16)
    else
        error("$(fmt.compression_code) is an unsupported compression code!")
    end
    if format == "double"
        samples = convert_to_double(samples)
    end
    samples
end

function write_pcm_samples(io::IO, fmt::WAVFormat, samples::AbstractArray{T}) where T <: Integer
    nbits = bits_per_sample(fmt)
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

function wavread(io::IO; subrange=Void, format="double")
    chunk_size = read_header(io)
    samples = Array{Float64, 1}()
    nbits = 0
    sample_rate = Float32(0.0)
    opt = Dict{Symbol, Any}()

    # Note: This assumes that the format chunk is written in the file before the data chunk. The
    # specification does not require this assumption, but most real files are written that way.

    # Subtract the size of the format field from chunk_size; now it holds the size
    # of all the sub-chunks
    chunk_size -= 4
    # GitHub Issue #18: Check if there is enough data to read another chunk
    subchunk_header_size = 4 + sizeof(UInt32)
    fmt = WAVFormat()
    while chunk_size >= subchunk_header_size
        # Read subchunk ID and size
        subchunk_id = Array{UInt8}(4)
        read!(io, subchunk_id)
        subchunk_size = read_le(io, UInt32)
        if subchunk_size > chunk_size
            chunk_size = 0
            break
        end
        chunk_size -= subchunk_header_size + subchunk_size
        # check the subchunk ID
        if subchunk_id == b"fmt "
            fmt = read_format(io, subchunk_size)
            sample_rate = Float32(fmt.sample_rate)
            nbits = bits_per_sample(fmt)
            opt[:fmt] = fmt
        elseif subchunk_id == b"data"
            if format == "size"
                return convert(Int, subchunk_size / fmt.block_align), convert(Int, fmt.nchannels)
            end
            samples = read_data(io, subchunk_size, fmt, format, make_range(subrange))
        else
            opt[Symbol(subchunk_id)] = Array{UInt8}(subchunk_size)
            read!(io, opt[Symbol(subchunk_id)])
        end
    end
    return samples, sample_rate, nbits, opt
end

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

function wavappend(samples::AbstractArray, io::IO)
    seekstart(io)
    chunk_size = read_header(io)
    subchunk_id = Array{UInt8}(4)
    read!(io, subchunk_id)
    subchunk_size = read_le(io, UInt32)
    if subchunk_id != b"fmt "
        error("First chunk is not the format")
    end
    fmt = read_format(io, subchunk_size)

    if fmt.nchannels != size(samples,2)
        error("Number of channels do not match")
    end

    data_length = size(samples, 1) * fmt.block_align

    seek(io,4)
    write_le(io, convert(UInt32, chunk_size + data_length))

    seek(io,64)
    subchunk_size = read_le(io, UInt32)
    seek(io,64)
    write_le(io, convert(UInt32, subchunk_size + data_length))

    seekend(io)
    write_data(io, fmt, samples)
end

function wavappend(samples::AbstractArray, filename::AbstractString)
    open(filename, "a+") do io
        wavappend(samples,io)
    end
end

wavwrite(y::AbstractArray, f::Real, filename::AbstractString) = wavwrite(y, filename, Fs=f)
wavwrite(y::AbstractArray, f::Real, n::Real, filename::AbstractString) = wavwrite(y, filename, Fs=f, nbits=n)

# support for writing native arrays...
wavwrite(y::AbstractArray{T}, io::IO) where T <: Integer = wavwrite(y, io, nbits=sizeof(T)*8)
wavwrite(y::AbstractArray{T}, filename::AbstractString) where T <: Integer = wavwrite(y, filename, nbits=sizeof(T)*8)
wavwrite(y::AbstractArray{Int32}, io::IO) = wavwrite(y, io, nbits=24)
wavwrite(y::AbstractArray{Int32}, filename::AbstractString) = wavwrite(y, filename, nbits=24)
wavwrite(y::AbstractArray{T}, io::IO) where T <: AbstractFloat = wavwrite(y, io, nbits=sizeof(T)*8, compression=WAVE_FORMAT_IEEE_FLOAT)
wavwrite(y::AbstractArray{T}, filename::AbstractString) where T <: AbstractFloat = wavwrite(y, filename, nbits=sizeof(T)*8, compression=WAVE_FORMAT_IEEE_FLOAT)

# FileIO integration support
load(s::Stream{format"WAV"}; kwargs...) = wavread(s.io; kwargs...)
save(s::Stream{format"WAV"}, data; kwargs...) = wavwrite(data, s.io; kwargs...)

load(f::File{format"WAV"}; kwargs...) = wavread(f.filename; kwargs...)
save(f::File{format"WAV"}, data; kwargs...) = wavwrite(data, f.filename; kwargs...)

end # module
