@testset "Writing WAV Files" begin
    @testset "Default WAVSink Properties" begin
        io = IOBuffer()
        WAVSink(io) do sink
            @test eltype(sink) == Float32
            @test nchannels(sink) == 2
            @test samplerate(sink) == 48000
        end
    end

    @testset "Lossless Round-Trip" begin
        freqs = [480, 960]
        sig = sin.(2pi*freqs' .* (0:99)/48000) * 0.3
        for T in [
                Float64,
                Float32,
                Fixed{Int64, 63},
                Fixed{Int32, 31},
                Fixed{Int32, 23},
                Fixed{Int32, 19},
                Fixed{Int16, 15},
                Fixed{Int8, 7},
                Fixed{Int8, 5},
                ]
            @testset "$T" begin
                testdata = map(T, sig)
                io = IOBuffer()
                WAVSink(io, eltype=T) do sink
                    write(sink, testdata)
                    @test eltype(sink) == T
                end
                seek(io, 0)
                WAVSource(io) do src
                    @test eltype(src) == T
                    @test read(src) == testdata
                end
            end
        end
    end
end
# These float array comparison functions are from dists.jl
# function absdiff(current::AbstractArray{T}, target::AbstractArray{T}) where T <: Real
#     @assert all(size(current) == size(target))
#     maximum(abs.(current - target))
# end

# function reldiff(current::T, target::T) where T <: Real
#     abs.((current - target)/(bool(target) ? target : 1))
# end

# function reldiff(current::AbstractArray{T}, target::AbstractArray{T}) where T <: Real
#     @assert all(size(current) == size(target))
#     maximum([reldiff(current[i], target[i]) for i in 1:numel(target)])
# end

# ## example from README, modified to use an IO buffer
# let
#     x = [0:7999;]
#     y = sin.(2 * pi * x / 8000)
#     io = IOBuffer()
#     WAV.wavwrite(y, io, Fs=8000)
#     seek(io, 0)
#     y, Fs = WAV.wavread(io)
#     y = cos.(2 * pi * x / 8000)
#     WAV.wavappend(y, io)
#     seek(io, 0)
#     y, Fs = WAV.wavread(io)
# end

# let
#     x = [0:7999;]
#     y = sin.(2 * pi * x / 8000)
#     WAV.wavwrite(y, "example.wav", Fs=8000)
#     y, Fs = WAV.wavread("example.wav")
#     y = cos.(2 * pi * x / 8000)
#     WAV.wavappend(y, "example.wav")
#     y, Fs = WAV.wavread("example.wav")
#     @test length(y) == (2 * length(x))
#     rm("example.wav")
# end

# ## default arguments, GitHub Issue #10
# let
#     tmp=rand(Float32,(10,2))
#     io = IOBuffer()
#     WAV.wavwrite(tmp, io; nbits=32)
#     seek(io, 0)
#     y, fs, nbits, extra = WAV.wavread(io; format="native")
#     @test typeof(y) == Array{Float32, 2}
#     @test fs == 8000.0
#     @test nbits == 32
#     @test WAV.isformat(extra[:fmt], WAV.WAVE_FORMAT_IEEE_FLOAT)
# end

# ## malformed subchunk header, GitHub Issue #18
# let
#     # Create a malformed WAV file
#     samples = rand(Float32, (10, 1))
#     io = IOBuffer()

#     const compression = WAV.get_default_compression(samples)
#     nbits = WAV.get_default_precision(samples, compression)

#     const nchannels = size(samples, 2)
#     const sample_rate = 8000
#     nbits = ceil(Integer, nbits / 8) * 8
#     const block_align = nbits / 8 * nchannels
#     const bps = sample_rate * block_align
#     const data_length::UInt32 = size(samples, 1) * block_align

#     const fmt = WAV.WAVFormat(compression,
#                               nchannels,
#                               sample_rate,
#                               bps,
#                               block_align,
#                               nbits,
#                               WAV.WAVFormatExtension())

#     WAV.write_header(io, UInt32(data_length + 37)) # 37 instead of 36 is the broken part
#     WAV.write_format(io, fmt)

#     # write the data subchunk header
#     WAV.write(io, b"data")
#     WAV.write_le(io, data_length) # UInt32
#     WAV.write_data(io, fmt, samples)

#     seek(io, 0)
#     y, fs, nbits, opt = WAV.wavread(io, format="native")
#     @test fs == 8000
#     @test nbits == 32
#     @test opt[:fmt].compression_code == compression
#     @test opt[:fmt].nchannels == nchannels
#     @test opt[:fmt].sample_rate == sample_rate
#     @test opt[:fmt].bytes_per_second == bps
#     @test opt[:fmt].block_align == block_align
#     @test WAV.isformat(opt[:fmt], compression)
#     @test WAV.bits_per_sample(opt[:fmt]) == nbits
#     @test samples == y
# end

# let
#     tmp=rand(Float64,(10,2))
#     io = IOBuffer()
#     WAV.wavwrite(tmp, io; nbits=64)
#     seek(io, 0)
#     y, fs, nbits, extra = WAV.wavread(io; format="native")
#     @test typeof(y) == Array{Float64, 2}
#     @test fs == 8000.0
#     @test nbits == 64
#     @test WAV.isformat(extra[:fmt], WAV.WAVE_FORMAT_IEEE_FLOAT)
#     @test WAV.isformat(extra[:fmt], WAV.WAVE_FORMAT_EXTENSIBLE)
# end

# let
#     tmp=rand(Float32,(10,2))
#     io = IOBuffer()
#     WAV.wavwrite(tmp, io; compression=WAV.WAVE_FORMAT_IEEE_FLOAT)
#     seek(io, 0)
#     y, fs, nbits, extra = WAV.wavread(io; format="native")
#     @test typeof(y) == Array{Float32, 2}
#     @test fs == 8000.0
#     @test nbits == 32
#     @test WAV.isformat(extra[:fmt], WAV.WAVE_FORMAT_IEEE_FLOAT)
# end

# ## Test wavread and wavwrite
# ## Generate some wav files for writing and reading
# for fs = (8000,11025,22050,44100,48000,96000,192000), nbits = (1,7,8,9,12,16,20,24,32,64), nsamples = convert(Array{Int}, [0; logspace(1, 4, 4)]), nchans = 1:4
#     ## Test wav files
#     ## The tolerance is based on the number of bits used to encode the file in wavwrite
#     tol = 2.0 / (2.0^(nbits - 1))

#     in_data = rand(nsamples, nchans)
#     if nsamples > 0
#         @test maximum(in_data) <= 1.0
#         @test minimum(in_data) >= -1.0
#     end
#     io = IOBuffer()
#     WAV.wavwrite(in_data, io, Fs=fs, nbits=nbits, compression=WAV.WAVE_FORMAT_PCM)
#     file_size = position(io)

#     ## Check for the common header identifiers
#     seek(io, 0)
#     @test read(io, UInt8, 4) == b"RIFF"
#     @test WAV.read_le(io, UInt32) == file_size - 8
#     @test read(io, UInt8, 4) == b"WAVE"

#     ## Check that wavread works on the wavwrite produced memory
#     seek(io, 0)
#     sz = WAV.wavread(io, format="size")
#     @test sz == (nsamples, nchans)

#     seek(io, 0)
#     out_data, out_fs, out_nbits, out_extra = WAV.wavread(io)
#     @test length(out_data) == nsamples * nchans
#     @test size(out_data, 1) == nsamples
#     @test size(out_data, 2) == nchans
#     @test typeof(out_data) == Array{Float64, 2}
#     @test out_fs == fs
#     @test out_nbits == nbits
#     fmt = out_extra[:fmt]
#     @test WAV.bits_per_sample(fmt) == nbits
#     @test WAV.isformat(fmt, WAV.WAVE_FORMAT_PCM)
#     @test fmt.nchannels == nchans
#     if nsamples > 0
#         @test absdiff(out_data, in_data) < tol
#     end
#     if nchans > 2
#         @test WAV.isformat(fmt, WAV.WAVE_FORMAT_EXTENSIBLE)
#     end

#     ## test the "subrange" option.
#     if nsamples > 0
#         seek(io, 0)
#         # Don't convert to Int, test if passing a float (nsamples/2) behaves as expected
#         subsamples = min(10, trunc(Int, nsamples / 2))
#         out_data, out_fs, out_nbits, out_extra = WAV.wavread(io, subrange=subsamples)
#         @test length(out_data) == subsamples * nchans
#         @test size(out_data, 1) == subsamples
#         @test size(out_data, 2) == nchans
#         @test typeof(out_data) == Array{Float64, 2}
#         @test out_fs == fs
#         @test out_nbits == nbits
#         fmt = out_extra[:fmt]
#         @test WAV.bits_per_sample(fmt) == nbits
#         @test WAV.isformat(fmt, WAV.WAVE_FORMAT_PCM)
#         @test fmt.nchannels == nchans
#         @test absdiff(out_data, in_data[1:subsamples, :]) < tol

#         seek(io, 0)
#         sr = convert(Int, min(5, trunc(Int, nsamples / 2))):convert(Int, min(23, nsamples - 1))
#         out_data, out_fs, out_nbits, out_extra = WAV.wavread(io, subrange=sr)
#         @test length(out_data) == length(sr) * nchans
#         @test size(out_data, 1) == length(sr)
#         @test size(out_data, 2) == nchans
#         @test typeof(out_data) == Array{Float64, 2}
#         @test out_fs == fs
#         @test out_nbits == nbits
#         fmt = out_extra[:fmt]
#         @test WAV.bits_per_sample(fmt) == nbits
#         @test WAV.isformat(fmt, WAV.WAVE_FORMAT_PCM)
#         @test fmt.nchannels == nchans
#         @test absdiff(out_data, in_data[sr, :]) < tol
#     end
# end

# ## Test wavappend
# x0=rand(2,3)
# x1=rand(3,3)
# x2=rand(4,3)
# io = IOBuffer()
# WAV.wavwrite(x0, io)
# WAV.wavappend(x1, io)
# WAV.wavappend(x2, io)
# seek(io, 0)
# x, fs, nbits, extra = WAV.wavread(io)
# @test x == [x0; x1; x2]

# ## Test native encoding of 8 bits
# for nchans = (1,2,4)
#     in_data_8 = reshape(typemin(UInt8):typemax(UInt8), (trunc(Int, 256 / nchans), nchans))
#     io = IOBuffer()
#     WAV.wavwrite(in_data_8, io)

#     seek(io, 0)
#     out_data_8, fs, nbits, extra = WAV.wavread(io, format="native")
#     @test fs == 8000
#     @test nbits == 8
#     fmt = extra[:fmt]
#     @test WAV.bits_per_sample(fmt) == nbits
#     @test WAV.isformat(fmt, WAV.WAVE_FORMAT_PCM)
#     @test fmt.nchannels == nchans
#     @test in_data_8 == out_data_8
# end

# ## Test native encoding of 16 bits
# for nchans = (1,2,4)
#     in_data_16 = reshape(typemin(Int16):typemax(Int16), (trunc(Int, 65536 / nchans), nchans))
#     io = IOBuffer()
#     WAV.wavwrite(in_data_16, io)

#     seek(io, 0)
#     out_data_16, fs, nbits, extra = WAV.wavread(io, format="native")
#     @test fs == 8000
#     @test nbits == 16
#     fmt = extra[:fmt]
#     @test WAV.bits_per_sample(fmt) == nbits
#     @test WAV.isformat(fmt, WAV.WAVE_FORMAT_PCM)
#     @test fmt.nchannels == nchans
#     @test in_data_16 == out_data_16
# end

# ## Test native encoding of 24 bits
# for nchans = (1,2,4)
#     in_data_24 = convert(Array{Int32}, reshape(-63:64, trunc(Int, 128 / nchans), nchans))
#     io = IOBuffer()
#     WAV.wavwrite(in_data_24, io)

#     seek(io, 0)
#     out_data_24, fs, nbits, extra = WAV.wavread(io, format="native")
#     @test fs == 8000
#     @test nbits == 24
#     fmt = extra[:fmt]
#     @test WAV.bits_per_sample(fmt) == nbits
#     @test WAV.isformat(fmt, WAV.WAVE_FORMAT_PCM)
#     @test fmt.nchannels == nchans
#     @test in_data_24 == out_data_24
# end

# ## Test encoding 32 bit values
# for nchans = (1,2,4)
#     in_data_single = convert(Array{Float32}, reshape(linspace(-1.0, 1.0, 128), trunc(Int, 128 / nchans), nchans))
#     io = IOBuffer()
#     WAV.wavwrite(in_data_single, io)

#     seek(io, 0)
#     out_data_single, fs, nbits, extra = WAV.wavread(io, format="native")
#     @test fs == 8000
#     @test nbits == 32
#     fmt = extra[:fmt]
#     @test WAV.bits_per_sample(fmt) == nbits
#     @test WAV.isformat(fmt, WAV.WAVE_FORMAT_IEEE_FLOAT)
#     @test fmt.nchannels == nchans
#     @test in_data_single == out_data_single
# end

# ## Test encoding 32 bit values outside the [-1, 1] range
# for nchans = (1,2,4)
#     nsamps = trunc(Int, 128 / nchans)
#     in_data_single = convert(Array{Float32}, reshape(-63:64, nsamps, nchans))
#     io = IOBuffer()
#     WAV.wavwrite(in_data_single, io)

#     seek(io, 0)
#     out_data_single, fs, nbits, extra = WAV.wavread(io, format="native")
#     @test fs == 8000
#     @test nbits == 32
#     fmt = extra[:fmt]
#     @test WAV.bits_per_sample(fmt) == nbits
#     @test WAV.isformat(fmt, WAV.WAVE_FORMAT_IEEE_FLOAT)
#     @test fmt.nchannels == nchans
#     @test [in_data_single[i, j] for i = 1:nsamps, j = 1:nchans] == out_data_single
# end

# ## Test encoding 64 bit values
# for nchans = (1,2,4)
#     in_data_single = convert(Array{Float64}, reshape(linspace(-1.0, 1.0, 128), trunc(Int, 128 / nchans), nchans))
#     io = IOBuffer()
#     WAV.wavwrite(in_data_single, io)

#     seek(io, 0)
#     out_data_single, fs, nbits, extra = WAV.wavread(io, format="native")
#     @test fs == 8000
#     @test nbits == 64
#     fmt = extra[:fmt]
#     @test WAV.bits_per_sample(fmt) == nbits
#     @test fmt.nchannels == nchans
#     @test WAV.isformat(fmt, WAV.WAVE_FORMAT_IEEE_FLOAT)
#     @test in_data_single == out_data_single
# end

# ## Test encoding 64 bit values outside the [-1, 1] range
# for nchans = (1,2,4)
#     nsamps = trunc(Int, 128 / nchans)
#     in_data_single = convert(Array{Float64}, reshape(-63:64, nsamps, nchans))
#     io = IOBuffer()
#     WAV.wavwrite(in_data_single, io)

#     seek(io, 0)
#     out_data_single, fs, nbits, extra = WAV.wavread(io, format="native")
#     @test fs == 8000
#     @test nbits == 64
#     fmt = extra[:fmt]
#     @test WAV.bits_per_sample(fmt) == nbits
#     @test WAV.isformat(fmt, WAV.WAVE_FORMAT_IEEE_FLOAT)
#     @test fmt.nchannels == nchans
#     @test [in_data_single[i, j] for i = 1:nsamps, j = 1:nchans] == out_data_single
# end

# ### Test A-Law and Mu-Law
# for nbits = (8, 16), nsamples = convert(Array{Int}, [0; logspace(1, 4, 4)]), nchans = 1:4, fmt=(WAV.WAVE_FORMAT_ALAW, WAV.WAVE_FORMAT_MULAW)
#     const fs = 8000.0
#     const tol = 2.0 / (2.0^6)
#     in_data = rand(nsamples, nchans)
#     if nsamples > 0
#         @test maximum(in_data) <= 1.0
#         @test minimum(in_data) >= -1.0
#     end
#     io = IOBuffer()
#     WAV.wavwrite(in_data, io, Fs=fs, nbits=nbits, compression=fmt)
#     file_size = position(io)

#     ## Check for the common header identifiers
#     seek(io, 0)
#     @test read(io, UInt8, 4) == b"RIFF"
#     @test WAV.read_le(io, UInt32) == file_size - 8
#     @test read(io, UInt8, 4) == b"WAVE"

#     ## Check that wavread works on the wavwrite produced memory
#     seek(io, 0)
#     sz = WAV.wavread(io, format="size")
#     @test sz == (nsamples, nchans)

#     seek(io, 0)
#     out_data, out_fs, out_nbits, out_extra = WAV.wavread(io)
#     @test length(out_data) == nsamples * nchans
#     @test size(out_data, 1) == nsamples
#     @test size(out_data, 2) == nchans
#     @test typeof(out_data) == Array{Float64, 2}
#     @test out_fs == fs
#     @test out_nbits == 8
#     @test WAV.bits_per_sample(out_extra[:fmt]) == 8
#     @test out_extra[:fmt].nchannels == nchans
#     @test WAV.isformat(out_extra[:fmt], fmt)
#     if nchans > 2
#         @test WAV.isformat(out_extra[:fmt], WAV.WAVE_FORMAT_EXTENSIBLE)
#     else
#         @test !WAV.isformat(out_extra[:fmt], WAV.WAVE_FORMAT_EXTENSIBLE)
#     end
#     if nsamples > 0
#         @test absdiff(out_data, in_data) < tol
#     end

#     ## test the "subrange" option.
#     if nsamples > 0
#         seek(io, 0)
#         # Don't convert to Int, test if passing a float (nsamples/2) behaves as expected
#         subsamples = min(10, trunc(Int, nsamples / 2))
#         out_data, out_fs, out_nbits, out_extra = WAV.wavread(io, subrange=subsamples)
#         @test length(out_data) == subsamples * nchans
#         @test size(out_data, 1) == subsamples
#         @test size(out_data, 2) == nchans
#         @test typeof(out_data) == Array{Float64, 2}
#         @test out_fs == fs
#         @test out_nbits == 8
#         @test WAV.bits_per_sample(out_extra[:fmt]) == 8
#         @test out_extra[:fmt].nchannels == nchans
#         @test WAV.isformat(out_extra[:fmt], fmt)
#         @test absdiff(out_data, in_data[1:subsamples, :]) < tol

#         seek(io, 0)
#         sr = convert(Int, min(5, trunc(Int, nsamples / 2))):convert(Int, min(23, nsamples - 1))
#         out_data, out_fs, out_nbits, out_extra = WAV.wavread(io, subrange=sr)
#         @test length(out_data) == length(sr) * nchans
#         @test size(out_data, 1) == length(sr)
#         @test size(out_data, 2) == nchans
#         @test typeof(out_data) == Array{Float64, 2}
#         @test out_fs == fs
#         @test out_nbits == 8
#         @test WAV.bits_per_sample(out_extra[:fmt]) == 8
#         @test out_extra[:fmt].nchannels == nchans
#         @test WAV.isformat(out_extra[:fmt], fmt)
#         @test absdiff(out_data, in_data[sr, :]) < tol
#     end
# end

# ### Test float formatting
# for nbits = (32, 64), nsamples = convert(Array{Int}, [0; logspace(1, 4, 4)]), nchans = 1:2, fmt=(WAV.WAVE_FORMAT_IEEE_FLOAT)
#     const fs = 8000.0
#     const tol = 1e-6
#     in_data = rand(nsamples, nchans)
#     if nsamples > 0
#         @test maximum(in_data) <= 1.0
#         @test minimum(in_data) >= -1.0
#     end
#     io = IOBuffer()
#     WAV.wavwrite(in_data, io, Fs=fs, nbits=nbits, compression=fmt)
#     file_size = position(io)

#     ## Check for the common header identifiers
#     seek(io, 0)
#     @test read(io, UInt8, 4) == b"RIFF"
#     @test WAV.read_le(io, UInt32) == file_size - 8
#     @test read(io, UInt8, 4) == b"WAVE"

#     ## Check that wavread works on the wavwrite produced memory
#     seek(io, 0)
#     sz = WAV.wavread(io, format="size")
#     @test sz == (nsamples, nchans)

#     seek(io, 0)
#     out_data, out_fs, out_nbits, out_extra = WAV.wavread(io)
#     @test length(out_data) == nsamples * nchans
#     @test size(out_data, 1) == nsamples
#     @test size(out_data, 2) == nchans
#     @test typeof(out_data) == Array{Float64, 2}
#     @test out_fs == fs
#     @test out_nbits == nbits
#     @test WAV.bits_per_sample(out_extra[:fmt]) == nbits
#     @test WAV.isformat(out_extra[:fmt], WAV.WAVE_FORMAT_IEEE_FLOAT)
#     @test out_extra[:fmt].nchannels == nchans
#     if nsamples > 0
#         @test absdiff(out_data, in_data) < tol
#     end

#     ## test the "subrange" option.
#     if nsamples > 0
#         seek(io, 0)
#         # Don't convert to Int, test if passing a float (nsamples/2) behaves as expected
#         subsamples = min(10, trunc(Int, nsamples / 2))
#         out_data, out_fs, out_nbits, out_extra = WAV.wavread(io, subrange=subsamples)
#         @test length(out_data) == subsamples * nchans
#         @test size(out_data, 1) == subsamples
#         @test size(out_data, 2) == nchans
#         @test typeof(out_data) == Array{Float64, 2}
#         @test out_fs == fs
#         @test out_nbits == nbits
#         @test WAV.bits_per_sample(out_extra[:fmt]) == nbits
#         @test out_extra[:fmt].nchannels == nchans
#         @test WAV.isformat(out_extra[:fmt], WAV.WAVE_FORMAT_IEEE_FLOAT)
#         @test absdiff(out_data, in_data[1:subsamples, :]) < tol

#         seek(io, 0)
#         sr = convert(Int, min(5, trunc(Int, nsamples / 2))):convert(Int, min(23, nsamples - 1))
#         out_data, out_fs, out_nbits, out_extra = WAV.wavread(io, subrange=sr)
#         @test length(out_data) == length(sr) * nchans
#         @test size(out_data, 1) == length(sr)
#         @test size(out_data, 2) == nchans
#         @test typeof(out_data) == Array{Float64, 2}
#         @test out_fs == fs
#         @test out_nbits == nbits
#         @test WAV.bits_per_sample(out_extra[:fmt]) == nbits
#         @test WAV.isformat(out_extra[:fmt], WAV.WAVE_FORMAT_IEEE_FLOAT)
#         @test out_extra[:fmt].nchannels == nchans
#         @test absdiff(out_data, in_data[sr, :]) < tol
#     end
# end

# ### Read unknown chunks
# let
#     const fs = 8000.0
#     in_data = rand(1024, 2)
#     io = IOBuffer()
#     in_chunks = Dict(:test=>[0x1, 0x2, 0x3])
#     WAV.wavwrite(in_data, io, Fs=fs, chunks=in_chunks)

#     seek(io, 0)
#     data, fs, nbits, ext = WAV.wavread(io)

#     @test haskey(ext, :test) == true
#     @test ext[:test] == in_chunks[:test]
# end

# ### playback
# # The playback tests don't work on travis.
# let
#     const fs = 44100.0
#     t = 1:44100;
#     in_data = sin.(5.0 * t / fs) * 1e-6;
#     #WAV.wavplay(in_data, fs);
#     #WAV.wavplay([in_data in_data], fs);
# end
