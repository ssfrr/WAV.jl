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
include("wavsource.jl")
include("wavsink.jl")

"""
    pushmetadata(opt, key, val)

Push the given metadata chunk `opt`, which is a dictionary of vectors of
metadata chunks.
"""
function pushmetadata(opt, key, val)
    symbkey = Symbol(strip(String(key)))
    if symbkey in keys(opt)
        push!(opt[symbkey], val)
    else
        opt[symbkey] = [val]
    end
end


end # module
