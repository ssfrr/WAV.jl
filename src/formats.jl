# Required WAV Chunk; The format chunk describes how the waveform data is stored
struct WAVFormat
    compression_code::UInt16
    nchannels::UInt16
    sample_rate::UInt32
    bytes_per_second::UInt32 # average bytes per second
    block_align::UInt16
    nbits::UInt16
    ext::WAVFormatExtension
end

# used by WAVE_FORMAT_EXTENSIBLE
struct WAVFormatExtension
    nbits::UInt16 # overrides nbits in WAVFormat type
    channel_mask::UInt32
    sub_format::Array{UInt8, 1} # 16 byte GUID

    WAVFormatExtension() = new(0, 0, Array{UInt8, 1}(0))
    WAVFormatExtension(nb, cm, sb) = new(nb, cm, sb)
end

function WAVFormatExtension(bytes)
    if isempty(bytes)
        return WAVFormatExtension()
    end
    # split bytes into valid_bits_per_sample, channel_mask, and sub_format
    valid_bits_per_sample = (convert(UInt16, bytes[2]) << 8) | convert(UInt16, bytes[1])
    channel_mask = (convert(UInt32, bytes[6]) << 24) | (convert(UInt32, bytes[5]) << 16) | (convert(UInt32, bytes[4]) << 8) | convert(UInt32, bytes[3])
    sub_format = bytes[7:end]
    return WAVFormatExtension(valid_bits_per_sample, channel_mask, sub_format)
end

const WAVE_FORMAT_PCM        = 0x0001 # PCM
const WAVE_FORMAT_IEEE_FLOAT = 0x0003 # IEEE float
const WAVE_FORMAT_ALAW       = 0x0006 # A-Law
const WAVE_FORMAT_MULAW      = 0x0007 # Mu-Law
const WAVE_FORMAT_EXTENSIBLE = 0xfffe # Extension!

isextensible(fmt::WAVFormat) = (fmt.compression_code == WAVE_FORMAT_EXTENSIBLE)
bits_per_sample(fmt::WAVFormat) = isextensible(fmt) ? fmt.ext.nbits : fmt.nbits

function read_format(io::IO, chunk_size::UInt32)
    # can I read in all of the fields at once?
    orig_chunk_size = convert(Int, chunk_size)
    if chunk_size < 16
        error("The WAVE Format chunk must be at least 16 bytes")
    end
    const compression_code = read_le(io, UInt16)
    const nchannels = read_le(io, UInt16)
    const sample_rate = read_le(io, UInt32)
    const bytes_per_second = read_le(io, UInt32)
    const block_align = read_le(io, UInt16)
    const nbits = read_le(io, UInt16)
    ext = Array{UInt8, 1}(0)
    chunk_size -= 16
    if chunk_size > 0
        const extra_bytes_length = read_le(io, UInt16)
        if extra_bytes_length == 22
            ext = read(io, UInt8, extra_bytes_length)
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


# DEFINE_GUIDSTRUCT("00000001-0000-0010-8000-00aa00389b71", KSDATAFORMAT_SUBTYPE_PCM);
const KSDATAFORMAT_SUBTYPE_PCM = [
0x01, 0x00, 0x00, 0x00,
0x00, 0x00,
0x10, 0x00,
0x80, 0x00,
0x00, 0xaa, 0x00, 0x38, 0x9b, 0x71
                                  ]
# DEFINE_GUIDSTRUCT("00000003-0000-0010-8000-00aa00389b71", KSDATAFORMAT_SUBTYPE_IEEE_FLOAT);
const KSDATAFORMAT_SUBTYPE_IEEE_FLOAT = [
0x03, 0x00, 0x00, 0x00,
0x00, 0x00,
0x10, 0x00,
0x80, 0x00,
0x00, 0xaa, 0x00, 0x38, 0x9b, 0x71
                                         ]
# DEFINE_GUID(KSDATAFORMAT_SUBTYPE_MULAW, 0x00000007, 0x0000, 0x0010, 0x80, 0x00, 0x00, 0xaa, 0x00, 0x38, 0x9b, 0x71);
const KSDATAFORMAT_SUBTYPE_MULAW = [
0x07, 0x00, 0x00, 0x00,
0x00, 0x00,
0x10, 0x00,
0x80, 0x00,
0x00, 0xaa, 0x00, 0x38, 0x9b, 0x71
                                    ]
#DEFINE_GUID(KSDATAFORMAT_SUBTYPE_ALAW, 0x00000006, 0x0000, 0x0010, 0x80, 0x00, 0x00, 0xaa, 0x00, 0x38, 0x9b, 0x71);
const KSDATAFORMAT_SUBTYPE_ALAW = [
0x06, 0x00, 0x00, 0x00,
0x00, 0x00,
0x10, 0x00,
0x80, 0x00,
0x00, 0xaa, 0x00, 0x38, 0x9b, 0x71
                                   ]

function isformat(fmt::WAVFormat, code)
    if code != WAVE_FORMAT_EXTENSIBLE && isextensible(fmt)
        subtype = Array{UInt8, 1}(0)
        if code == WAVE_FORMAT_PCM
            subtype = KSDATAFORMAT_SUBTYPE_PCM
        elseif code == WAVE_FORMAT_IEEE_FLOAT
            subtype = KSDATAFORMAT_SUBTYPE_IEEE_FLOAT
        elseif code == WAVE_FORMAT_ALAW
            subtype = KSDATAFORMAT_SUBTYPE_ALAW
        elseif code == WAVE_FORMAT_MULAW
            subtype = KSDATAFORMAT_SUBTYPE_MULAW
        else
            return false
        end
        return subtype == fmt.ext.sub_format
    end
    return fmt.compression_code == code
end
