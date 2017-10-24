# used by WAVE_FORMAT_EXTENSIBLE
struct WAVFormatExtension
    nbits::UInt16 # overrides nbits in WAVFormat type
    channel_mask::UInt32
    sub_format::Array{UInt8, 1} # 16 byte GUID

    WAVFormatExtension() = new(0, 0, Vector{UInt8}(0))
    WAVFormatExtension(nb, cm, sb) = new(nb, cm, sb)
end

# read a WAVFormatExtension from raw bytes
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

const COMPRESSIONS = (:none, :mulaw, :alaw)

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

function WAVFormat(eltype, nchannels, samplerate, compression)
    fmtcode, ext_fmtcode = formatcodes(eltype, nchannels, compression)
    containerbytes = ceil(Int, validbits(eltype) / 8)
    blockalign = containerbytes * nchannels
    bps = samplerate * blockalign
    ext = if needs_extensible(eltype, nchannels)
            WAVFormatExtension(validbits(eltype), 0, ext_fmtcode)
        else
            WAVFormatExtension()
        end

    WAVFormat(fmtcode, nchannels, samplerate, bps, blockalign,
              containerbytes*8, ext)
end

function Base.show(io::IO, fmt::WAVFormat)
    println(io, "WAVFormat")
    println(io, "  compression_code: $(formatname(fmt.compression_code))")
    println(io, "  nchannels: $(fmt.nchannels)")
    println(io, "  sample_rate: $(fmt.sample_rate)")
    println(io, "  bits_per_second: $(fmt.bytes_per_second)")
    println(io, "  block_align: $(fmt.block_align)")
    println(io, "  nbits: $(fmt.nbits)")
    if isextensible(fmt)
        println(io, "  extended data:")
        println(io, "    valid_bits_per_sample: $(fmt.ext.valid_bits_per_sample)")
        println(io, "    channel_mask: $(fmt.ext.channel_mask)")
        println(io, "    sub_format: $(subformatname(fmt.ext.valid_bits_per_sample))")
    end
end

function needs_extensible(eltype, nchannels)
    # criteria from http://www-mmsp.ece.mcgill.ca/Documents/AudioFormats/WAVE/WAVE.html
    if nchannels > 2
        true
    elseif eltype <: Fixed && (validbits(eltype) > 16 || validbits(eltype) % 8 != 0)
        true
    else
        false
    end
end

"""
    formatcodes(eltype, channels, compression)

Returns the normal and extended format codes appropriate for the given
parameters.
"""
function formatcodes(eltype, nchannels, compression)
    if needs_extensible(eltype, nchannels)
        if eltype <: AbstractFloat && compression == :none
            WAVE_FORMAT_EXTENSIBLE, KSDATAFORMAT_SUBTYPE_IEEE_FLOAT
        elseif eltype <: Fixed && compression == :none
            WAVE_FORMAT_EXTENSIBLE, KSDATAFORMAT_SUBTYPE_PCM
        elseif eltype == Fixed{Int16, 15} && compression == :alaw
            WAVE_FORMAT_EXTENSIBLE, KSDATAFORMAT_SUBTYPE_ALAW
        elseif eltype == Fixed{Int16, 15} && compression == :mulaw
            WAVE_FORMAT_EXTENSIBLE, KSDATAFORMAT_SUBTYPE_MULAW
        else
            error("Invalid eltype, nchannel, compression combination: $eltype, $nchannels, $compression")
        end
    else
        if eltype <: AbstractFloat && compression == :none
            WAVE_FORMAT_IEEE_FLOAT, nothing
        elseif eltype <: Fixed && compression == :none
            WAVE_FORMAT_PCM, nothing
        elseif eltype == Fixed{Int16, 15} && compression == :alaw
            WAVE_FORMAT_ALAW, nothing
        elseif eltype == Fixed{Int16, 15} && compression == :mulaw
            WAVE_FORMAT_MULAW, nothing
        else
            error("Invalid eltype, nchannel, compression combination: $eltype, $nchannels, $compression")
        end
    end
end

const WAVE_FORMAT_PCM        = 0x0001 # PCM
const WAVE_FORMAT_IEEE_FLOAT = 0x0003 # IEEE float
const WAVE_FORMAT_ALAW       = 0x0006 # A-Law
const WAVE_FORMAT_MULAW      = 0x0007 # Mu-Law
const WAVE_FORMAT_EXTENSIBLE = 0xfffe # Extension!

function formatname(code)
    if code == WAVE_FORMAT_PCM
        "WAVE_FORMAT_PCM"
    elseif code == WAVE_FORMAT_IEEE_FLOAT
        "WAVE_FORMAT_IEEE_FLOAT"
    elseif code == WAVE_FORMAT_IEEE_ALAW
        "WAVE_FORMAT_IEEE_ALAW"
    elseif code == WAVE_FORMAT_IEEE_MULAW
        "WAVE_FORMAT_IEEE_MULAW"
    elseif code == WAVE_FORMAT_IEEE_EXTENSIBLE
        "WAVE_FORMAT_IEEE_EXTENSIBLE"
    else
        "Unknown format $code"
    end
end

function subformatname(code)
    if code == KSDATAFORMAT_SUBTYPE_PCM
        "KSDATAFORMAT_SUBTYPE_PCM"
    elseif code == KSDATAFORMAT_SUBTYPE_IEEE_FLOAT
        "KSDATAFORMAT_SUBTYPE_IEEE_FLOAT"
    elseif code == KSDATAFORMAT_SUBTYPE_IEEE_ALAW
        "KSDATAFORMAT_SUBTYPE_IEEE_ALAW"
    elseif code == KSDATAFORMAT_SUBTYPE_IEEE_MULAW
        "KSDATAFORMAT_SUBTYPE_IEEE_MULAW"
    elseif code == KSDATAFORMAT_SUBTYPE_IEEE_EXTENSIBLE
        "KSDATAFORMAT_SUBTYPE_IEEE_EXTENSIBLE"
    else
        "Unknown subformat $code"
    end
end

isextensible(fmt::WAVFormat) = (fmt.compression_code == WAVE_FORMAT_EXTENSIBLE)

"""
    validbits(x::WAVFormat)
    validbits(x::Type)

Return the number of valid bits per sample in the given WAV format block or
Julia type.
"""
validbits(fmt::WAVFormat) = isextensible(fmt) ? fmt.ext.nbits : fmt.nbits
validbits(::Type{<:Fixed{T, f}}) where {T, f} = f+1
validbits(::Type{Float32}) = 32
validbits(::Type{Float64}) = 64

function read_format(io::IO, chunkbytes, bytesleft)
    if chunkbytes < 16
        error("The WAVE Format chunk must be at least 16 bytes")
    end
    compression_code = read_le(io, UInt16)
    nchannels = read_le(io, UInt16)
    sample_rate = read_le(io, UInt32)
    bytes_per_second = read_le(io, UInt32)
    block_align = read_le(io, UInt16)
    nbits = read_le(io, UInt16)
    ext = Vector{UInt8}(0)
    bytesleft -= 16
    chunkbytes -= 16
    if chunkbytes > 0
        extra_bytes_length = read_le(io, UInt16)
        bytesleft -= 2
        if extra_bytes_length == 22
            ext = read(io, UInt8, extra_bytes_length)
            length(ext) == extra_bytes_length || throw(EOFError())
            bytesleft -= extra_bytes_length
        end
    end
    return WAVFormat(compression_code,
                     nchannels,
                     sample_rate,
                     bytes_per_second,
                     block_align,
                     nbits,
                     WAVFormatExtension(ext)),
           bytesleft
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

    # return the number of bytes written
    len+8
end


# DEFINE_GUIDSTRUCT("00000001-0000-0010-8000-00aa00389b71", KSDATAFORMAT_SUBTYPE_PCM);
const KSDATAFORMAT_SUBTYPE_PCM = [
    0x01, 0x00, 0x00, 0x00,
    0x00, 0x00,
    0x10, 0x00,
    0x80, 0x00,
    0x00, 0xaa, 0x00, 0x38, 0x9b, 0x71]
# DEFINE_GUIDSTRUCT("00000003-0000-0010-8000-00aa00389b71", KSDATAFORMAT_SUBTYPE_IEEE_FLOAT);
const KSDATAFORMAT_SUBTYPE_IEEE_FLOAT = [
    0x03, 0x00, 0x00, 0x00,
    0x00, 0x00,
    0x10, 0x00,
    0x80, 0x00,
    0x00, 0xaa, 0x00, 0x38, 0x9b, 0x71]
# DEFINE_GUID(KSDATAFORMAT_SUBTYPE_MULAW, 0x00000007, 0x0000, 0x0010, 0x80, 0x00, 0x00, 0xaa, 0x00, 0x38, 0x9b, 0x71);
const KSDATAFORMAT_SUBTYPE_MULAW = [
    0x07, 0x00, 0x00, 0x00,
    0x00, 0x00,
    0x10, 0x00,
    0x80, 0x00,
    0x00, 0xaa, 0x00, 0x38, 0x9b, 0x71]
#DEFINE_GUID(KSDATAFORMAT_SUBTYPE_ALAW, 0x00000006, 0x0000, 0x0010, 0x80, 0x00, 0x00, 0xaa, 0x00, 0x38, 0x9b, 0x71);
const KSDATAFORMAT_SUBTYPE_ALAW = [
    0x06, 0x00, 0x00, 0x00,
    0x00, 0x00,
    0x10, 0x00,
    0x80, 0x00,
    0x00, 0xaa, 0x00, 0x38, 0x9b, 0x71]

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
