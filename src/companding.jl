decodemulaw(encoded) = reinterpret(Fixed{Int16, 15}, mudecode_table[encoded+1])
decodealaw(encoded) = reinterpret(Fixed{Int16, 15}, adecode_table[encoded+1])

function encodemulaw(sample)
    const cBias = 0x84
    const cClip = 32635

    const sampleSign = (sample >>> 8) & 0x80
    if sampleSign != 0
        sample = -sample
    end
    if sample > cClip
        sample = cClip
    end
    sample = sample + cBias
    const sampleExponent = muencode_table[(sample >>> 7) + 1]
    const mantissa = (sample >> (sampleExponent+3)) & 0x0F
    UInt8((~ (sampleSign | (sampleExponent << 4) | mantissa)) & 0xff)
end

function encodealaw(sample)
    const cBias = 0x84
    const cClip = 32635
    const sampleSign = ((~sample >>> 8) & 0x80)
    if sampleSign == 0
        sample = -sample
    end
    if sample > cClip
        sample = cClip
    end
    compressedByte = 0
    if sample >= 256
        const sampleExponent = aencode_table[((sample >>> 8) & 0x7f) + 1]
        const mantissa = (sample >>> (sampleExponent + 3) ) & 0x0f
        compressedByte = ((sampleExponent << 4) | mantissa) & 0xff
    else
        compressedByte = (sample >>> 4) & 0xff
    end
    compressedByte = xor(xor(sampleSign, 0x55), compressedByte)
    UInt8(compressedByte & 0xff)
end

# Quantized μ-law algorithm -- Use a look up table to convert
# From Wikipedia, ITU-T Recommendation G.711 and G.191 specify the following intervals:
#
# ---------------------------------------+--------------------------------
#  14 bit Binary Linear input code       | 8 bit Compressed code
# ---------------------------------------+--------------------------------
# +8158 to +4063 in 16 intervals of 256  |  0x80 + interval number
# +4062 to +2015 in 16 intervals of 128  |  0x90 + interval number
# +2014 to +991 in 16 intervals of 64    |  0xA0 + interval number
# +990 to +479 in 16 intervals of 32     |  0xB0 + interval number
# +478 to +223 in 16 intervals of 16     |  0xC0 + interval number
# +222 to +95 in 16 intervals of 8       |  0xD0 + interval number
# +94 to +31 in 16 intervals of 4        |  0xE0 + interval number
# +30 to +1 in 15 intervals of 2         |  0xF0 + interval number
# 0                                      |  0xFF
# −1                                     |  0x7F
# −31 to −2 in 15 intervals of 2         |  0x70 + interval number
# −95 to −32 in 16 intervals of 4        |  0x60 + interval number
# −223 to −96 in 16 intervals of 8       |  0x50 + interval number
# −479 to −224 in 16 intervals of 16     |  0x40 + interval number
# −991 to −480 in 16 intervals of 32     |  0x30 + interval number
# −2015 to −992 in 16 intervals of 64    |  0x20 + interval number
# −4063 to −2016 in 16 intervals of 128  |  0x10 + interval number
# −8159 to −4064 in 16 intervals of 256  |  0x00 + interval number
# ---------------------------------------+--------------------------------
const mudecode_table = Int16[
    -32124,-31100,-30076,-29052,-28028,-27004,-25980,-24956,
    -23932,-22908,-21884,-20860,-19836,-18812,-17788,-16764,
    -15996,-15484,-14972,-14460,-13948,-13436,-12924,-12412,
    -11900,-11388,-10876,-10364, -9852, -9340, -8828, -8316,
    -7932, -7676, -7420, -7164, -6908, -6652, -6396, -6140,
    -5884, -5628, -5372, -5116, -4860, -4604, -4348, -4092,
    -3900, -3772, -3644, -3516, -3388, -3260, -3132, -3004,
    -2876, -2748, -2620, -2492, -2364, -2236, -2108, -1980,
    -1884, -1820, -1756, -1692, -1628, -1564, -1500, -1436,
    -1372, -1308, -1244, -1180, -1116, -1052,  -988,  -924,
    -876,  -844,  -812,  -780,  -748,  -716,  -684,  -652,
    -620,  -588,  -556,  -524,  -492,  -460,  -428,  -396,
    -372,  -356,  -340,  -324,  -308,  -292,  -276,  -260,
    -244,  -228,  -212,  -196,  -180,  -164,  -148,  -132,
    -120,  -112,  -104,   -96,   -88,   -80,   -72,   -64,
    -56,   -48,   -40,   -32,   -24,   -16,    -8,     -1,
    32124, 31100, 30076, 29052, 28028, 27004, 25980, 24956,
    23932, 22908, 21884, 20860, 19836, 18812, 17788, 16764,
    15996, 15484, 14972, 14460, 13948, 13436, 12924, 12412,
    11900, 11388, 10876, 10364,  9852,  9340,  8828,  8316,
    7932,  7676,  7420,  7164,  6908,  6652,  6396,  6140,
    5884,  5628,  5372,  5116,  4860,  4604,  4348,  4092,
    3900,  3772,  3644,  3516,  3388,  3260,  3132,  3004,
    2876,  2748,  2620,  2492,  2364,  2236,  2108,  1980,
    1884,  1820,  1756,  1692,  1628,  1564,  1500,  1436,
    1372,  1308,  1244,  1180,  1116,  1052,   988,   924,
    876,   844,   812,   780,   748,   716,   684,   652,
    620,   588,   556,   524,   492,   460,   428,   396,
    372,   356,   340,   324,   308,   292,   276,   260,
    244,   228,   212,   196,   180,   164,   148,   132,
    120,   112,   104,    96,    88,    80,    72,    64,
    56,    48,    40,    32,    24,    16,     8,     0]

@assert length(mudecode_table) == 256

const muencode_table = [
    0,0,1,1,2,2,2,2,3,3,3,3,3,3,3,3,
    4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
    5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
    5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
    6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
    6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
    6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
    6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
    7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,
    7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,
    7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,
    7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,
    7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,
    7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,
    7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,
    7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7]

@assert length(muencode_table) == 256

# Quantized A-law algorithm -- Use a look up table to convert
const adecode_table = Int16[
    -5504, -5248, -6016, -5760, -4480, -4224, -4992, -4736,
    -7552, -7296, -8064, -7808, -6528, -6272, -7040, -6784,
    -2752, -2624, -3008, -2880, -2240, -2112, -2496, -2368,
    -3776, -3648, -4032, -3904, -3264, -3136, -3520, -3392,
    -22016,-20992,-24064,-23040,-17920,-16896,-19968,-18944,
    -30208,-29184,-32256,-31232,-26112,-25088,-28160,-27136,
    -11008,-10496,-12032,-11520,-8960, -8448, -9984, -9472,
    -15104,-14592,-16128,-15616,-13056,-12544,-14080,-13568,
    -344,  -328,  -376,  -360,  -280,  -264,  -312,  -296,
    -472,  -456,  -504,  -488,  -408,  -392,  -440,  -424,
    -88,   -72,   -120,  -104,  -24,   -8,    -56,   -40,
    -216,  -200,  -248,  -232,  -152,  -136,  -184,  -168,
    -1376, -1312, -1504, -1440, -1120, -1056, -1248, -1184,
    -1888, -1824, -2016, -1952, -1632, -1568, -1760, -1696,
    -688,  -656,  -752,  -720,  -560,  -528,  -624,  -592,
    -944,  -912,  -1008, -976,  -816,  -784,  -880,  -848,
    5504,  5248,  6016,  5760,  4480,  4224,  4992,  4736,
    7552,  7296,  8064,  7808,  6528,  6272,  7040,  6784,
    2752,  2624,  3008,  2880,  2240,  2112,  2496,  2368,
    3776,  3648,  4032,  3904,  3264,  3136,  3520,  3392,
    22016, 20992, 24064, 23040, 17920, 16896, 19968, 18944,
    30208, 29184, 32256, 31232, 26112, 25088, 28160, 27136,
    11008, 10496, 12032, 11520, 8960,  8448,  9984,  9472,
    15104, 14592, 16128, 15616, 13056, 12544, 14080, 13568,
    344,   328,   376,   360,   280,   264,   312,   296,
    472,   456,   504,   488,   408,   392,   440,   424,
    88,    72,   120,   104,    24,     8,    56,    40,
    216,   200,   248,   232,   152,   136,   184,   168,
    1376,  1312,  1504,  1440,  1120,  1056,  1248,  1184,
    1888,  1824,  2016,  1952,  1632,  1568,  1760,  1696,
    688,   656,   752,   720,   560,   528,   624,   592,
    944,   912,  1008,   976,   816,   784,   880,   848]

@assert length(adecode_table) == 256

const aencode_table = [
    1,1,2,2,3,3,3,3,
    4,4,4,4,4,4,4,4,
    5,5,5,5,5,5,5,5,
    5,5,5,5,5,5,5,5,
    6,6,6,6,6,6,6,6,
    6,6,6,6,6,6,6,6,
    6,6,6,6,6,6,6,6,
    6,6,6,6,6,6,6,6,
    7,7,7,7,7,7,7,7,
    7,7,7,7,7,7,7,7,
    7,7,7,7,7,7,7,7,
    7,7,7,7,7,7,7,7,
    7,7,7,7,7,7,7,7,
    7,7,7,7,7,7,7,7,
    7,7,7,7,7,7,7,7,
    7,7,7,7,7,7,7,7]

@assert length(aencode_table) == 128