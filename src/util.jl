"""
    pushmetadata(opt, key, val)

Push the given metadata chunk `opt`, which is a dictionary of vectors of
metadata chunks, Keyed on the 4-byte chunk ID, as a `Vector{UInt8}`, e.g.
`b"data"`
"""
function pushmetadata(opt, key::Vector{UInt8}, val)
    if key in keys(opt)
        push!(opt[key], val)
    else
        opt[key] = [val]
    end
end

"""
Convert the underlying wav format to a Julia type that we can
use to parameterize sources/sinks for type-stable reads and writes.

Companded streams are represented as 16-bit fixed-point because that's the
user-facing representation.
"""
function type_from_format(fmt)
    nbits = validbits(fmt)
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
            error("FORMAT_IEEE_FLOAT only supports 32 and 64-bit floats (not $nbits)")
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
