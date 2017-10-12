# used by WAVE_FORMAT_EXTENSIBLE
struct WAVFormatExtension
    nbits::UInt16 # overrides nbits in WAVFormat type
    channel_mask::UInt32
    sub_format::Array{UInt8, 1} # 16 byte GUID
    
    WAVFormatExtension() = new(0, 0, Array{UInt8, 1}(0))
    WAVFormatExtension(nb, cm, sb) = new(nb, cm, sb)
end

# Required WAV Chunk; The format chunk describes how the waveform data is stored
struct WAVFormat
    compression_code::UInt16
    nchannels::UInt16
    sample_rate::UInt32
    bytes_per_second::UInt32 # average bytes per second
    block_align::UInt16
    nbits::UInt16
    ext::WAVFormatExtension

    WAVFormat() = new(WAVE_FORMAT_PCM, 0, 0, 0, 16, 0, WAVFormatExtension())
    WAVFormat(cc, nchan, fs, bps, ba, nb, e) = new(cc, nchan, fs, bps, ba, nb, e)
end

const WAVE_FORMAT_PCM        = 0x0001 # PCM
const WAVE_FORMAT_IEEE_FLOAT = 0x0003 # IEEE float
const WAVE_FORMAT_ALAW       = 0x0006 # A-Law
const WAVE_FORMAT_MULAW      = 0x0007 # Mu-Law
const WAVE_FORMAT_EXTENSIBLE = 0xfffe # Extension!

isextensible(fmt::WAVFormat) = (fmt.compression_code == WAVE_FORMAT_EXTENSIBLE)
bits_per_sample(fmt::WAVFormat) = isextensible(fmt) ? fmt.ext.nbits : fmt.nbits

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
