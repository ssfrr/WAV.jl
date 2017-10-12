WAV.jl
======

[![WAV](http://pkg.julialang.org/badges/WAV_0.4.svg)](http://pkg.julialang.org/?pkg=WAV&ver=0.4)
[![WAV](http://pkg.julialang.org/badges/WAV_0.5.svg)](http://pkg.julialang.org/?pkg=WAV&ver=0.5)
[![WAV](http://pkg.julialang.org/badges/WAV_0.6.svg)](http://pkg.julialang.org/?pkg=WAV&ver=0.6)
[![Build Status](https://travis-ci.org/dancasimiro/WAV.jl.png)](https://travis-ci.org/dancasimiro/WAV.jl)
[![Windows Build Status](https://ci.appveyor.com/api/projects/status/github/dancasimiro/wav.jl?branch=master&svg=true)](https://ci.appveyor.com/project/dancasimiro/wav-jl)
[![Coverage Status](https://coveralls.io/repos/dancasimiro/WAV.jl/badge.png)](https://coveralls.io/r/dancasimiro/WAV.jl)

This is a pure Julia package to read and write the WAV audio file format.

Getting Started
---------------

WAV provides an interface to work with `.wav` files using the [FileIO](https://github.com/JuliaIO/FileIO.jl) interface, and `SampleBuf` objects from the [SampledSignals](https://github.com/JuliaAudio/SampledSignals.jl) package, implemented in pure Julia. Here is an example to get you started. It
generates some data, writes it to a file and then reads the data back. Data is kept in its native format (sample rate and type) when possible, but automatically converted by the SampledSignals infrastructure when necessary. Additionally, if the WAV file is 8-bit with μ-law or a-law companding, it will be automatically decoded to 16-bit integer values.

```jlcon
julia> using WAV
julia> x = [0:7999;]
julia> y = sin.(2 * pi * x / 8000)
julia> save("example.wav", y, samplerate=8000)
julia> y = load("example.wav")
julia> y *= 2
```

For MATLAB compatibility we also provide a `wavread` function with the same arguments.

load
----

This function reads the samples from a WAV file, returning a `SampleBuf`.

```julia
function wavread(io::IO; subrange=Any, format="double")
function wavread(filename::String; subrange=Any, format="double")
```

The available options, and the default values, are:

The returned values are:

* ``y``: The acoustic samples; A matrix is returned for files that contain multiple channels.
* ``Fs``: The sampling frequency
* ``nbits``: The number of bits used to encode each sample
* ``opt``: A ```Dict{Symbol, Any}``` of optional chunks found in the WAV file.

The dictionary returned in the ``opt`` field depends on the contents
of the WAV file. All valid WAV files will contain a "fmt" chunk. The
"fmt" entry in the dictionary will contain an instance of type
``WAVFormat``. The ``WAVFormat`` type is defined as::

```julia
immutable WAVFormat
    compression_code::UInt16
    nchannels::UInt16
    sample_rate::UInt32
    bytes_per_second::UInt32 # average bytes per second
    block_align::UInt16
    nbits::UInt16
    ext::WAVFormatExtension
end
```

The ```ext``` field of type ```WAVFormatExtension``` is defined as::

```julia
immutable WAVFormatExtension
    nbits::UInt16 # overrides nbits in WAVFormat type
    channel_mask::UInt32
    sub_format::Array{UInt8, 1} # 16 byte GUID
    WAVFormatExtension() = new(0, 0, Array(UInt8, 0))
    WAVFormatExtension(nb, cm, sb) = new(nb, cm, sb)
end
```

You can use the ```isformat``` function to test how the samples are
encoded, without worrying about the ```WAVFormatExtension```
type. Extended WAV files were added to deal with some ambiguity in the
original specification.

```julia
isformat(fmt::WAVFormat, code)
```

The ```isformat``` function takes the format object from the ```opt```
output dictionary of ```wavread``` and one of ```WAV_FORMAT_```
constants, respectively. The function returns ```true``` when the
samples are encoded in the specified ```code```.

The following functions are also defined to make this function compatible with MATLAB:

```julia
wavread(filename::String, fmt::String) = wavread(filename, format=fmt)
wavread(filename::String, N::Int) = wavread(filename, subrange=N)
wavread(filename::String, N::Range1{Int}) = wavread(filename, subrange=N)
wavread(filename::String, N::Int, fmt::String) = wavread(filename, subrange=N, format=fmt)
wavread(filename::String, N::Range1{Int}, fmt::String) = wavread(filename, subrange=N, format=fmt)
```

wavwrite
--------

Writes samples to a RIFF/WAVE file io object. The ``io`` argument
accepts either an ``IO`` object or a filename (``String``). The
function assumes that the sample rate is 8 kHz and uses 16 bits to
encode each sample. Both of these values can be changed with the
options parameter. Each column of the data represents a different
channel. Stereo files should contain two columns. The options are
passed via an ``Options`` object (see the :ref:`options page
<options-module>`).

```julia
function wavwrite(samples::Array, io::IO; Fs=8000, nbits=16, compression=WAVE_FORMAT_PCM)
function wavwrite(samples::Array, filename::String; Fs=8000, nbits=16, compression=WAVE_FORMAT_PCM)
```

The available options, and the default values, are:

   * ``Fs`` (default = ``8000``): sampling frequency
   * ``nbits`` (default = ``16``): number of bits used to encode each
     sample
   * ``compression (default = ``WAV_FORMAT_PCM``)``: controls the type of encoding used in the file

The type of the input array, samples, also affects the generated
file. "Native" WAVE files are written when integers are passed into
wavwrite. This means that the literal values are written into the
file. The input ranges are as follows for integer samples.

| N Bits | y Data Type | y Data Range           | Output Format |
|--------|-------------|------------------------|---------------|
| 8      | uint8       | 0 <= y <= 255          | uint8         |
| 16     | int16       | –32768 <= y <= +32767  | int16         |
| 24     | int32       | –2^23 <= y <= 2^23 – 1 | int32         |

If samples contains floating point values, the input data ranges
are the following.

| N Bits | y Data Type      | y Data Range       | Output Format |
|--------|------------------|--------------------|---------------|
| 8      | single or double |  –1.0 <= y < +1.0  | uint8         |
| 16     | single or double |  –1.0 <= y < +1.0  | int16         |
| 24     | single or double |  –1.0 <= y < +1.0  | int32         |
| 32     | single or double |  –1.0 <= y <= +1.0 | single        |

Floating point (single and double precision) values are written to the
file unaltered. The library will not modify the data range or representation.

The following functions are also defined to make this function
compatible with MATLAB:

```julia
wavwrite(y::Array, f::Real, filename::String) = wavwrite(y, filename, Fs=f)
wavwrite(y::Array, f::Real, N::Real, filename::String) = wavwrite(y, filename, Fs=f, nbits=N)
wavwrite{T<:Integer}(y::Array{T}, io::IO) = wavwrite(y, io, nbits=sizeof(T)*8)
wavwrite{T<:Integer}(y::Array{T}, filename::String) = wavwrite(y, filename, nbits=sizeof(T)*8)
wavwrite(y::Array{Int32}, io::IO) = wavwrite(y, io, nbits=24)
wavwrite(y::Array{Int32}, filename::String) = wavwrite(y, filename, nbits=24)
wavwrite{T<:FloatingPoint}(y::Array{T}, io::IO) = wavwrite(y, io, nbits=sizeof(T)*8, compression=WAVE_FORMAT_IEEE_FLOAT)
wavwrite{T<:FloatingPoint}(y::Array{T}, filename::String) = wavwrite(y, filename, nbits=sizeof(T)*8, compression=WAVE_FORMAT_IEEE_FLOAT)
```

wavappend
---------

Append samples to an existing WAV file.  All parameters (data type and range, output format, number of bits, number of channels, etc.) are assumed to match.

```julia
function wavappend(samples::Array, io::IO)
function wavappend(samples::Array, filename::String)
```

Related Packages
----------------

[FLAC.jl](https://github.com/dmbates/FLAC.jl) includes
an ```mmap``` based WAV [reader](https://github.com/dmbates/FLAC.jl/blob/master/src/WAV.jl).
