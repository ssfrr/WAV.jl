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