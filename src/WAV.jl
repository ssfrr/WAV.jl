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

# chunks are keyed by their ID, e.g. `b"data"`. We do this rather than using
# `Symbol`s so we can capture spaces, e.g. `b"cue "`
const ChunkKey = Vector{UInt8}
const SUBCHUNK_HEADER_SIZE = 8

include("util.jl")
include("formats.jl")
include("companding.jl")
include("fileio.jl")
include("wavsource.jl")
include("wavsink.jl")

# implement do syntax to close stream automatically
for StrType in [:WAVSource, :WAVSink]
    @eval function $StrType(fn::Function, io::IO; kwargs...)
        stream = $StrType(io; kwargs...)
        try
            fn(stream)
        finally
            close(stream)
        end
    end
end

end # module
