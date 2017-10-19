The files in this directory are for testing WAV.jl's WAV reading capability.

## `48000_stereo_X.wav`

These files were created in Audacity. They are each 100-sample stereo files at 48000Hz sampling rate. The left channel is a 480Hz sin wav, the right is 960Hz, both with an amplitude of 0.3. They were created as native float32 samples and then exported as each of the given sample types.

## `48000_stereo_pcm16_long_filelength.wav`

This file is the same as the first set of files, except that the length field of the `RIFF` header has been modified to be longer than it actually is (and in fact extends past the end of the file)

## `48000_stereo_pcm16_long_datalength.wav`

This file is the same as the first set of files, except that the length field of the `data` chunk has been modified to be longer than it actually is (and in fact extends past the end of the file)

## `Utopia Critical Stop.WAV`

This file came from http://www-mmsp.ece.mcgill.ca/Documents/AudioFormats/WAVE/Samples.html, and has a `fact` chunk that comes after the `data` chunk. The fact chunk contains the four characters `FILT` and not the number of samples. It also contains several other post-data chunks, including two `DISP` chunks and a `LIST` chunk with an `INFO` chunk with `ICMT`, `ICOP`, and `ISBJ` subchunks. It is from the system directory `C:\WINNT\Media` on a Windows 2000 system.

## `Ptjunk.wav`

This file came from http://www-mmsp.ece.mcgill.ca/Documents/AudioFormats/WAVE/Samples.html, and has trailing garbage after the `data` chunk. It is 9 mu-law-encoded samples at 8000Hz with an odd length intermediate chunk (type XxXx, with a padding byte). Note it does not have a `fact` chunk, which is required for any non-PCM files.

## `Pmiscck.wav`
This file came from http://www-mmsp.ece.mcgill.ca/Documents/AudioFormats/WAVE/Samples.html, it is 9 mu-law-encoded samples at 8000Hz with an odd length intermediate chunk (type XxXx, with a padding byte). Note it does not have a `fact` chunk, which is required for any non-PCM files.