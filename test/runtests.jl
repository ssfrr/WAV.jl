## -*-Julia-*-
## Test suite for Julia's WAV module
using WAV
using Base.Test: @testset, @test
using TestSetExtensions

import FileIO
using FileIO: load, save

using SampledSignals: SampleBuf, nchannels, samplerate, nframes, metadata
using SampledSignals: PCM8Sample, PCM16Sample, PCM20Sample, PCM24Sample, PCM32Sample, PCM64Sample

using Suppressor: @capture_err, @color_output

@testset ExtendedTestSet "WAV.jl Tests" begin
    include("reading.jl")
    include("writing.jl")
end