@testset "File Reading" begin
    AUDIODIR = joinpath(@__DIR__, "audio")
    expected_frames = 100
    expected_channels = 2
    expected_samplerate = 48000

    # re-create the signal in the wav files. They were originally generated
    # in audacity, which natively uses 32-bit floats
    freqs = [480, 960]
    expected_data = map(Float32, sin.(2pi*freqs' .* (0:99)/48000) * 0.3)

    for (T, filename) in [
            (Float32, "48000_stereo_float32.wav"),
            (Float64, "48000_stereo_float64.wav"),
            (PCM8Sample, "48000_stereo_pcm8.wav"),
            (PCM16Sample, "48000_stereo_pcm16.wav"),
            (PCM20Sample, "48000_stereo_pcm20.wav"),
            (PCM24Sample, "48000_stereo_pcm24.wav"),
            (PCM32Sample, "48000_stereo_pcm32.wav"),
            (PCM16Sample, "48000_stereo_pcm16_long_filelength.wav"),
            (PCM16Sample, "48000_stereo_pcm16_long_datalength.wav"),
            ]
        loaded = WAV.load(FileIO.query(joinpath(AUDIODIR, filename)))
        @test loaded isa SampleBuf
        @test eltype(loaded) == T
        @test nchannels(loaded) == expected_channels
        @test samplerate(loaded) == expected_samplerate
        @test nframes(loaded) == expected_frames
        @test loaded ≈ map(T, expected_data)
    end

    # these ones need a somewhat lower precision then implied by their 16-bit
    # datatype because they're companded from 8-bit data
    for (T, filename) in [
            (PCM16Sample, "48000_stereo_alaw.wav"),
            (PCM16Sample, "48000_stereo_ulaw.wav"),
            ]
        loaded = WAV.load(FileIO.query(joinpath(AUDIODIR, filename)))
        @test loaded isa SampleBuf
        @test eltype(loaded) == T
        @test nchannels(loaded) == expected_channels
        @test samplerate(loaded) == expected_samplerate
        @test nframes(loaded) == expected_frames
        @test isapprox(loaded, map(T, expected_data), rtol=0.045)
    end
end