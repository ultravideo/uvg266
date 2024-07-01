uvg266
=======
An open-source VVC encoder licensed under 3-clause BSD

Join channel [#ultravideo](https://web.libera.chat/#ultravideo) in [Libera.Chat](https://libera.chat/) IRC network to contact us or come to our Discord [![Discord](https://img.shields.io/discord/973260924288901140?style=plastic)](https://discord.gg/fZpub7BPUA)

uvg266 is still under development. Speed and RD-quality will continue to improve.

https://ultravideo.fi/uvg266.html for more information.

- Linux [![uvg266_tests](https://github.com/ultravideo/uvg266/actions/workflows/uvg266.yml/badge.svg)](https://github.com/ultravideo/uvg266/actions/workflows/uvg266.yml)
- Windows [![Build status](https://ci.appveyor.com/api/projects/status/c1gwnnkyt5lycqka?svg=true)](https://ci.appveyor.com/project/Ultravideo/uvg266)

## Table of Contents

- [Using uvg266](#using-uvg266)
  - [Example:](#example)
  - [Parameters](#parameters)
  - [LP-GOP syntax](#lp-gop-syntax)
- [Presets](#presets)
- [uvg266 library](#uvg266-library)
- [Compiling uvg266](#compiling-uvg266)
  - [CMake](#cmake)
  - [Docker](#docker)
- [Paper](#paper)
- [Contributing to uvg266](#contributing-to-uvg266)
  - [Code documentation](#code-documentation)
  - [For version control we try to follow these conventions:](#for-version-control-we-try-to-follow-these-conventions)
  - [Testing](#testing)
  - [Unit tests](#unit-tests)
  - [Code style](#code-style)

## Using uvg266 VVC

### Debugging:

    ./uvg266 -i BQMall_832x480_60.yuv -o BQMall.266 -n 10 --no-sao --threads=0 --no-wpp -p 1 --rd=0 --fast-residual-cost=32 --no-deblock > debug.txt
    
    ./DecoderAnalyserApp -b BQMall.266 --TraceFile=trace.txt --TraceRule=D_COMMON,D_CABAC,D_SYNTAX,D_NALUNITHEADER,D_HEADER:poc>=0 -o rec.yuv

### Example:

    uvg266 --input BQMall_832x480_60.yuv --output out.vvc

The mandatory parameters are input and output. If the resolution of the input file is not in the filename, or when pipe is used, the input resolution must also be given: ```--input-res=1920x1080```.

The default input format is 8-bit yuv420p for 8-bit and yuv420p10le for 10-bit. Input format and bitdepth can be selected with ```--input-format``` and ```--input-bitdepth```.

Speed and compression quality can be selected with ```--preset```, or by setting the options manually.

### Parameters

[comment]: # (BEGIN UVG266 HELP MESSAGE)
```
uvg266 -i <input> --input-res <width>x<height> -o <output>

Required:
  -i, --input <filename>     : Input file
      --input-res <res>      : Input resolution [auto]
                                   - auto: Detect from file name.
                                   - <int>x<int>: width times height
  -o, --output <filename>    : Output file

Presets:
      --preset <preset>      : Set options to a preset [medium]
                                   - ultrafast, superfast, veryfast, faster,
                                     fast, medium, slow, slower, veryslow
                                     placebo

Input:
  -n, --frames <integer>     : Number of frames to code [all]
      --seek <integer>       : First frame to code [0]
      --input-fps <num>[/<denom>] : Frame rate of the input video [25]
      --source-scan-type <string> : Source scan type [progressive]
                                   - progressive: Progressive scan
                                   - tff: Top field first
                                   - bff: Bottom field first
      --input-format <string> : P420 or P400 [P420]
      --input-bitdepth <int> : 8-16 [8]
      --loop-input           : Re-read input file forever.
      --input-file-format <string> : Input file format [auto]
                                    - auto: Check the file ending for format
                                    - y4m (skips frame headers)
                                    - yuv

Options:
      --help                 : Print this help message and exit.
      --version              : Print version information and exit.
      --(no-)aud             : Use access unit delimiters. [disabled]
      --debug <filename>     : Output internal reconstruction.
      --(no-)cpuid           : Enable runtime CPU optimizations. [enabled]
      --hash <string>        : Decoded picture hash [checksum]
                                   - none: 0 bytes
                                   - checksum: 18 bytes
                                   - md5: 56 bytes
      --(no-)psnr            : Calculate PSNR for frames. [enabled]
      --(no-)info            : Add encoder info SEI. [enabled]
      --stats-file-prefix    : A prefix used for stats files that include
                               bits, lambda, distortion, and qp for each ctu.
                               These are meant for debugging and are not
                               written unless the prefix is defined.
      --cabac-debug-file     : A debug file for cabac context.
                               Ignore this, it is only for tests.

Video structure:
  -q, --qp <integer>         : Quantization parameter [22]
  -p, --period <integer>     : Period of intra pictures [64]
                                   - 0: Only first picture is intra.
                                   - 1: All pictures are intra.
                                   - N: Every Nth picture is intra.
      --vps-period <integer> : How often the video parameter set is re-sent [0]
                                   - 0: Only send VPS with the first frame.
                                   - N: Send VPS with every Nth intra frame.
  -r, --ref <integer>        : Number of reference frames, in range 1..15 [4]
      --gop <string>         : GOP structure [lp-g4d3t1]
                                   -  0: Disabled
                                   -  8: B-frame pyramid of length 8
                                   - 16: B-frame pyramid of length 16
                                   - lp-<string>: Low-delay P/B-frame GOP
                                     (e.g. lp-g8d4t2, see README)
      --intra-qp-offset <int>: QP offset for intra frames [-51..51] [auto]
                                   - N: Set QP offset to N.
                                   - auto: Select offset automatically based
                                     on GOP length.
      --(no-)open-gop        : Use open GOP configuration. [enabled]
      --cqmfile <filename>   : Read custom quantization matrices from a file.
      --scaling-list <string>: Set scaling list mode. [off]
                                   - off: Disable scaling lists.
                                   - custom: use custom list (with --cqmfile).
                                   - default: Use default lists.
      --bitrate <integer>    : Target bitrate [0]
                                   - 0: Disable rate control.
                                   - N: Target N bits per second.
      --rc-algorithm <string>: Select used rc-algorithm. [lambda]
                                   - lambda: rate control from:
                                     DOI: 10.1109/TIP.2014.2336550 
                                   - oba: DOI: 10.1109/TCSVT.2016.2589878
      --(no-)intra-bits      : Use Hadamard cost based allocation for intra
                               frames. Default on for gop 8 and off for lp-gop
      --(no-)clip-neighbour  : On oba based rate control whether to clip 
                               lambda values to same frame's ctus or previous'.
                               Default on for RA GOPS and disabled for LP.
      --(no-)lossless        : Use lossless coding. [disabled]
      --mv-constraint <string> : Constrain movement vectors. [none]
                                   - none: No constraint
                                   - frametile: Constrain within the tile.
                                   - frametilemargin: Constrain even more.
      --roi <filename>       : Use a delta QP map for region of interest.
                               Reads an array of delta QP values from a file.
                               Text and binary files are supported and detected
                               from the file extension (.txt/.bin). If a known
                               extension is not found, the file is treated as
                               a text file. The file can include one or many
                               ROI frames each in the following format:
                               width and height of the QP delta map followed
                               by width * height delta QP values in raster
                               order. In binary format, width and height are
                               32-bit integers whereas the delta QP values are
                               signed 8-bit values. The map can be of any size
                               and will be scaled to the video size. The file
                               reading will loop if end of the file is reached.
                               See roi.txt in the examples folder.
      --set-qp-in-cu         : Set QP at CU level keeping pic_init_qp_minus26.
                               in PPS and slice_qp_delta in slize header zero.
      --(no-)erp-aqp         : Use adaptive QP for 360 degree video with
                               equirectangular projection. [disabled]
      --level <number>       : Use the given HEVC level in the output and give
                               an error if level limits are exceeded. [6.2]
                                   - 1, 2, 2.1, 3, 3.1, 4, 4.1, 5, 5.1, 5.2, 6,
                                     6.1, 6.2
      --force-level <number> : Same as --level but warnings instead of errors.
      --high-tier            : Used with --level. Use high tier bitrate limits
                               instead of the main tier limits during encoding.
                               High tier requires level 4 or higher.
      --(no-)vaq <integer>   : Enable variance adaptive quantization with given
                               strength, in range 1..20. Recommended: 5.
                               [disabled]
      --chroma-qp-in         : List of input values used for mapping the luma
                               QP into chroma qp. [17,27,32,44]
      --chroma-qp-out        : List of output values used for mapping the luma
                               QP into chroma qp. These two lists have to be
                               same length, start with same value, and can
                               contain maximum 16 or 36 - starting value
                               elements. [17,27,32,44]
      --(no-)dual-tree       : Use separate CTU structure for luma and
                               chroma in intra slices.

Compression tools:
      --(no-)deblock <beta:tc> : Deblocking filter. [0:0]
                                   - beta: Between -6 and 6
                                   - tc: Between -6 and 6
      --sao <string>         : Sample Adaptive Offset [full]
                                   - off: SAO disabled
                                   - band: Band offset only
                                   - edge: Edge offset only
                                   - full: Full SAO
      --alf <string>         : Adaptive Loop Filter [off]
                                   - off: ALF disabled
                                   - no-cc: ALF enabled without cross component
                                            refinement
                                   - full: Full ALF
      --(no-)rdoq            : Rate-distortion optimized quantization [enabled]
      --(no-)rdoq-skip       : Skip RDOQ for 4x4 blocks. [disabled]
      --(no-)dep-quant       : Use dependent quantization. [disabled]
      --(no-)signhide        : Sign hiding [disabled]
      --rd <integer>         : Intra mode search complexity [0]
                                   - 0: Skip intra if inter is good enough.
                                   - 1: Rough intra mode search with SATD.
                                   - 2: Refine intra mode search with SSE.
                                   - 3: Enable intra chroma mode search.
                                   - 4: Try all intra modes.
      --(no-)mv-rdo          : Rate-distortion optimized motion vector costs
                               [disabled]
      --(no-)zero-coeff-rdo  : If a CU is set inter, check if forcing zero
                               residual improves the RD cost. [enabled]
      --(no-)full-intra-search : Try all intra modes during rough search.
                               [disabled]
      --(no-)transform-skip  : Try transform skip [disabled]
      --(no-)chroma-transform-skip : Try transform skip for chroma 
                                     blocks. [disabled]
      --tr-skip-max-size     : Max log2 size of transform skip 2..5 [2]
      --me <string>          : Integer motion estimation algorithm [hexbs]
                                   - hexbs: Hexagon Based Search
                                   - tz:    Test Zone Search
                                   - full:  Full Search
                                   - full8, full16, full32, full64
                                   - dia:   Diamond Search
      --me-steps <integer>   : Motion estimation search step limit. Only
                               affects 'hexbs' and 'dia'. [-1]
      --subme <integer>      : Fractional pixel motion estimation level [4]
                                   - 0: Integer motion estimation only
                                   - 1: + 1/2-pixel horizontal and vertical
                                   - 2: + 1/2-pixel diagonal
                                   - 3: + 1/4-pixel horizontal and vertical
                                   - 4: + 1/4-pixel diagonal
      --pu-depth-inter <int>-<int> : Maximum and minimum split depths where
                                     inter search is performed 0..8. [0-3]
                                   - Accepts a list of values separated by ','
                                     for setting separate depths per GOP layer
                                     (values can be omitted to use the first
                                     value for the respective layer).
      --pu-depth-intra <int>-<int> : Maximum and minimum split depths where
                                     intra search is performed 0..8. [1-4]
                                   - Accepts a list of values separated by ','
                                     for setting separate depths per GOP layer
                                     (values can be omitted to use the first
                                     value for the respective layer).
      --ml-pu-depth-intra    : Predict the pu-depth-intra using machine
                                learning trees, overrides the
                                --pu-depth-intra parameter. [disabled]
      --mtt-depth-intra      : Depth of mtt for intra slices 0..3.[0]
      --mtt-depth-intra-chroma : Depth of mtt for chroma dual tree in
                                      intra slices 0..3.[0]
      --mtt-depth-inter      : Depth of mtt for inter slices 0..3.[0]
                              All MTTs are currently experimental and
                              require disabling some avx2 optimizations.
      --max-bt-size          : maximum size for a CU resulting from
                                   a bt split. A singular value shared for all
                                   or a list of three values for the different
                                   slices types (intra, inter, intra-chroma)
                                   can be provided. [64, 64, 32]
      --max-tt-size          : maximum size for a CU resulting from
                                   a tt split. A singular value shared for all
                                   or a list of three values for the different
                                   slices types (intra, inter, intra-chroma)
                                   can be provided. [64, 64, 32]
      --intra-rough-granularity : How many levels are used for the
                                   logarithmic intra rough search. 0..4
                                   With 0 all of the modes are checked 
                                   in a single level, 1 checks every second
                                   mode is checked on first level and then
                                   second level checks the modes surrounding
                                   the three best modes. [2]
      --(no-)combine-intra-cus: Whether the encoder tries to code a cu
                                   on lower depth even when search is not
                                   performed on said depth. Should only
                                   be disabled if cus absolutely must not
                                   be larger than limited by the search.
                                   [enabled]
      --force-inter          : Force the encoder to use inter always.
                               This is mostly for debugging and is not
                               guaranteed to produce sensible bitstream or
                               work at all. [disabled]
      --(no-)bipred          : Bi-prediction [disabled]
      --cu-split-termination <string> : CU split search termination [zero]
                                   - off: Don't terminate early.
                                   - zero: Terminate when residual is zero.
      --me-early-termination <string> : Motion estimation termination [on]
                                   - off: Don't terminate early.
                                   - on: Terminate early.
                                   - sensitive: Terminate even earlier.
      --fast-residual-cost <int> : Skip CABAC cost for residual coefficients
                                   when QP is below the limit. [0]
      --fast-coeff-table <string> : Read custom weights for residual
                                    coefficients from a file instead of using
                                    defaults [default]
      --fast-rd-sampling : Enable learning data sampling for fast coefficient
                           table generation
      --fastrd-accuracy-check : Evaluate the accuracy of fast coefficient
                                prediction
      --fastrd-outdir : Directory to which to output sampled data or accuracy
                        data, into <fastrd-outdir>/0.txt to 50.txt, one file
                        for each QP that blocks were estimated on
      --(no-)intra-rdo-et    : Check intra modes in rdo stage only until
                               a zero coefficient CU is found. [disabled]
      --(no-)early-skip      : Try to find skip cu from merge candidates.
                               Perform no further search if skip is found.
                               For rd = 0..1: Try the first candidate.
                               For rd = 2.. : Try the best candidate based
                                              on luma satd cost. [enabled]
      --max-merge <integer>  : Maximum number of merge candidates, 1..6 [6]
      --(no-)implicit-rdpcm  : Implicit residual DPCM. Currently only supported
                               with lossless coding. [disabled]
      --(no-)tmvp            : Temporal motion vector prediction [enabled]
      --(no-)mrl             : Enable use of multiple reference lines in intra
                               predictions.
      --(no-)mip             : Enable matrix weighted intra prediction.
      --(no-)lfnst           : Enable low frequency non-separable transform.
                                 [disabled]
      --(no-)isp             : Enable intra sub partitions. [disabled]
                               Experimental, requires disabling some avx2
                               optimizations.
      --mts <string>         : Multiple Transform Selection [off].
                               (Currently only implemented for intra
                               and has effect only when rd >= 2)
                                   - off: MTS disabled
                                   - intra: MTS applied only for intra blocks.
                                   - inter: MTS applied only for inter blocks.
                                   - both: MTS applied for both intra and inter
                                           blocks.
                                   - implicit: uses implicit MTS. Applies DST7
                                               instead of DCT2 to certain intra
                                               blocks.
      --(no-)jccr            : Joint coding of chroma residual.
                               Requires rdo >= 2. [disabled]
      --(no-)cclm            : Cross component linear model. 
                               Extra chroma prediction modes that are formed
                               via linear transformation from the luma
                               prediction. Requires rdo >= 3. [disabled
      --(no-)amvr            : Adaptive Motion Vector Resolution.
                               Code some MVs with reduced resolution [disabled]

Parallel processing:
      --threads <integer>    : Number of threads to use [auto]
                                   - 0: Process everything with main thread.
                                   - N: Use N threads for encoding.
                                   - auto: Select automatically.
      --owf <integer>        : Frame-level parallelism [auto]
                                   - N: Process N+1 frames at a time.
                                   - auto: Select automatically.
      --(no-)wpp             : Wavefront parallel processing. [enabled]
                               Enabling tiles automatically disables WPP.
                               To enable WPP with tiles, re-enable it after
                               enabling tiles. Enabling wpp with tiles is,
                               however, an experimental feature since it is
                               not supported in any HEVC profile.
      --tiles <int>x<int>    : Split picture into width x height uniform tiles.
      --tiles-width-split <string>|u<int> :
                                   - <string>: A comma-separated list of tile
                                               column pixel coordinates.
                                   - u<int>: Number of tile columns of uniform
                                             width.
      --tiles-height-split <string>|u<int> :
                                   - <string>: A comma-separated list of tile
                                               row column pixel coordinates.
                                   - u<int>: Number of tile rows of uniform
                                             height.
      --slices <string>      : Control how slices are used.
                                   - tiles: Put tiles in independent slices.
                                   - wpp: Put rows in dependent slices.
                                   - tiles+wpp: Do both.
      --partial-coding <x-offset>!<y-offset>!<slice-width>!<slice-height>
                             : Encode partial frame.
                               Parts must be merged to form a valid bitstream.
                               X and Y are CTU offsets.
                               Slice width and height must be divisible by CTU
                               in pixels unless it is the last CTU row/column.
                               This parameter is used by kvaShare.

Video Usability Information:
      --sar <width:height>   : Specify sample aspect ratio
      --overscan <string>    : Specify crop overscan setting [undef]
                                   - undef, show, crop
      --videoformat <string> : Specify video format [undef]
                                   - undef, component, pal, ntsc, secam, mac
      --range <string>       : Specify color range [tv]
                                   - tv, pc
      --colorprim <string>   : Specify color primaries [undef]
                                   - undef, bt709, bt470m, bt470bg,
                                     smpte170m, smpte240m, film, bt2020
      --transfer <string>    : Specify transfer characteristics [undef]
                                   - undef, bt709, bt470m, bt470bg,
                                     smpte170m, smpte240m, linear, log100,
                                     log316, iec61966-2-4, bt1361e,
                                     iec61966-2-1, bt2020-10, bt2020-12
      --colormatrix <string> : Specify color matrix setting [undef]
                                   - undef, bt709, fcc, bt470bg, smpte170m,
                                     smpte240m, GBR, YCgCo, bt2020nc, bt2020c
      --chromaloc <integer>  : Specify chroma sample location (0 to 5) [0]
```
[comment]: # (END UVG266 HELP MESSAGE)


### LP-GOP syntax
The LP-GOP syntax is "lp-g(num)d(num)t(num)", where
- g = GOP length.
- d = Number of GOP layers.
- t = How many references to skip for temporal scaling, where 4 means only
  every fourth picture needs to be decoded.

```
QP
+4  o o o o  
+3   o   o    o o o o
+2     o       o o o    ooooooo
+1 o       o o       o o       o ooooooooo
   g8d4t1    g8d3t1    g8d2t1    g8d1t1
```

## Presets
The names of the presets are the same as with x264: ultrafast,
superfast, veryfast, faster, fast, medium, slow, slower, veryslow and
placebo. The effects of the presets are listed in the following table,
where the names have been abbreviated to fit the layout in GitHub.

|                      | 0-uf  | 1-sf  | 2-vf  | 3-fr  | 4-f   | 5-m   | 6-s   | 7-sr  | 8-vs  | 9-p   |
| -------------------- | ----- | ----- | ----- | ----- | ----- | ----- | ----- | ----- | ----- | ----- |
| rd                   | 0     | 0     | 0     | 0     | 0     | 0     | 1     | 2     | 2     | 2     |
| pu-depth-intra       | 2-3   | 2-3   | 2-3   | 2-3   | 1-3   | 1-4   | 1-4   | 1-4   | 1-4   | 1-4   |
| pu-depth-inter       | 1-2   | 1-2   | 1-3   | 1-3   | 1-3   | 0-3   | 0-3   | 0-3   | 0-3   | 0-3   |
| me                   | hexbs | hexbs | hexbs | hexbs | hexbs | hexbs | hexbs | hexbs | tz    | tz    |
| gop                  | 8     | 8     | 8     | 8     | 8     | 16    | 16    | 16    | 16    | 16    |
| ref                  | 1     | 1     | 1     | 1     | 2     | 4     | 4     | 4     | 4     | 4     |
| bipred               | 1     | 1     | 1     | 1     | 1     | 1     | 1     | 1     | 1     | 1     |
| deblock              | 1     | 1     | 1     | 1     | 1     | 1     | 1     | 1     | 1     | 1     |
| signhide             | 0     | 0     | 0     | 0     | 0     | 0     | 0     | 1     | 1     | 1     |
| subme                | 0     | 2     | 2     | 4     | 4     | 4     | 4     | 4     | 4     | 4     |
| sao                  | off   | full  | full  | full  | full  | full  | full  | full  | full  | full  |
| rdoq                 | 0     | 0     | 0     | 0     | 0     | 1     | 1     | 1     | 1     | 1     |
| rdoq-skip            | 0     | 0     | 0     | 0     | 0     | 0     | 0     | 0     | 0     | 0     |
| transform-skip       | 0     | 0     | 0     | 0     | 0     | 0     | 0     | 0     | 1     | 1     |
| mv-rdo               | 0     | 0     | 0     | 0     | 0     | 0     | 0     | 0     | 0     | 1     |
| full-intra-search    | 0     | 0     | 0     | 0     | 0     | 0     | 0     | 0     | 0     | 0     |
| cu-split-termination | zero  | zero  | zero  | zero  | zero  | zero  | zero  | zero  | zero  | off   |
| me-early-termination | sens. | sens. | sens. | sens. | sens. | on    | on    | off   | off   | off   |
| intra-rdo-et         | 0     | 0     | 0     | 0     | 0     | 0     | 0     | 0     | 0     | 0     |
| early-skip           | 1     | 1     | 1     | 1     | 1     | 1     | 1     | 1     | 1     | 0     |
| fast-residual-cost   | 28    | 28    | 28    | 0     | 0     | 0     | 0     | 0     | 0     | 0     |
| max-merge            | 5     | 5     | 5     | 5     | 5     | 5     | 5     | 5     | 5     | 5     |
| cclm                 | 0     | 0     | 0     | 0     | 0     | 0     | 0     | 0     | 1     | 1     |
| dual-tree            | 0     | 0     | 0     | 0     | 0     | 0     | 0     | 0     | 1     | 1     |
| jccr                 | 0     | 0     | 0     | 0     | 0     | 0     | 0     | 0     | 1     | 1     |
| mip                  | 0     | 0     | 0     | 0     | 0     | 0     | 0     | 0     | 1     | 1     |
| mrl                  | 0     | 0     | 0     | 0     | 0     | 0     | 0     | 0     | 1     | 1     |
| mts                  | off   | off   | off   | off   | off   | off   | off   | off   | both  | both  |
| dep-quant            | 0     | 0     | 0     | 0     | 0     | 0     | 0     | 0     | 1     | 1     |

## uvg266 library
See [uvg266.h](src/uvg266.h) for the library API and its
documentation.

When using the static uvg266 library on Windows, macro `UVG_STATIC_LIB`
must be defined. On other platforms it's not strictly required.

The needed linker and compiler flags can be obtained with pkg-config.


## Compiling uvg266
If you have trouble regarding compiling the source code, please make an
[issue](https://github.com/ultravideo/uvg266/issues) about in Github.
Others might encounter the same problem and there is probably much to
improve in the build process. We want to make this as simple as
possible.

**CMakeLists.txt assumes that the x86 based CPUs are 64bit and support AVX2**

### CMake
Depending on the platform, some additional tools are required for compiling uvg266 with CMake.
For Ubuntu, the required packages are `build-essential cmake`.

Run the following commands to generate the build scripts

    cd build
    cmake ..
    
Then depending on your generator settings you might want to use Make to compile and install uvg266, force with `-G 'Unix Makefiles'` in CMake command

    make
    sudo make install
    
or Ninja, force with `-G Ninja` in the CMake command

    ninja
    sudo ninja install

Visual Studio natively supports opening the `CMakeLists.txt` of the CMake build package has been installed.
Otherwise CMake-CLI can be used to generate the Visual Studio project files.
**When building shared library with visual studio/MSys2/MinGW the tests will fail to link, so they are disabled**

### Docker
This project includes a [Dockerfile](./Dockerfile), which enables building for Docker.
Build using Docker: `docker build -t uvg266 .`
Example usage: `docker run -i -a STDIN -a STDOUT uvg266 -i - --input-res=320x240 -o - < testfile_320x240.yuv > out.266`
For other examples, see [Dockerfile](./Dockerfile)

## Paper

**The paper is open access**

Please cite [this paper](https://ieeexplore.ieee.org/document/9690938) for uvg266:

```M. Viitanen, J. Sainio, A. Mercat, A. Lemmetti, and J. Vanne, “From HEVC to VVC: the First Development Steps of a Practical Intra Video Encoder,” Accepted to IEEE Transactions on Consumer Electronics```

Or in BibTex:

```
@ARTICLE{uvg266_2022, 
  author={Viitanen, Marko and Sainio, Joose and Mercat, Alexandre and Lemmetti, Ari and Vanne, Jarno},
  journal={IEEE Transactions on Consumer Electronics}, 
  title={From HEVC to VVC: the First Development Steps of a Practical Intra Video Encoder}, 
  year={2022},
  volume={},
  number={},
  doi={10.1109/TCE.2022.3146016}}
}
```

## Contributing to uvg266
We are happy to look at pull requests in Github. There is still lots of work to be done.


### Code documentation
You can generate Doxygen documentation pages by running the command
"doxygen docs.doxy". Here is a rough sketch of the module structure:
![uvg266 module hierarchy](https://github.com/ultravideo/uvg266/blob/master/doc/uvg266_module_hierarchy.png)


### For version control we try to follow these conventions:
- Master branch always produces a working bitstream (can be decoded with
  VTM).
- Commits for new features and major changes/fixes put to a sensibly
  named feature branch first and later merged to the master branch.
- Always merge the feature branch to the master branch, not the other
  way around, with fast-forwarding disabled if necessary.
- Every commit should at least compile.


### Testing
- Main automatic way of testing is with Github Actions. Commits, branches
  and pull requests are tested automatically.
  - Uninitialized variables and such are checked with Valgrind.
  - Bitstream validity is checked with VTM.
  - Compilation is checked on GCC and Clang on Linux, and Clang on OSX.
- Windows msys2 and msvc builds are checked automatically on Appveyor.
- If your changes change the bitstream, decode with VTM to check that
  it doesn't throw checksum errors or asserts.
- If your changes shouldn't alter the bitstream, check that they don't.
- Automatic compression quality testing is in the works.


### Unit tests
- There are some unit tests located in the tests directory. We would
  like to have more.
- The Visual Studio project links the unit tests against the actual .lib
  file used by the encoder. There is no Makefile as of yet.
- The unit tests use "greatest" unit testing framework. It is included
  as a submodule, but getting it requires the following commands to be
  run in the root directory of uvg266:

        git submodule init
        git submodule update

- On Linux, run ```make test```.


### Code style
We try to follow the following conventions:
- C99 without features not supported by Visual Studio 2015 (VLAs).
 - // comments allowed and encouraged.
- Follow overall conventions already established in the code.
- Indent by 2 spaces. (no tabs)
- { on the same line for control logic and on the next line for functions
- Reference and deference next to the variable name.
- Variable names in lowered characters with words divided by underscore.
- Maximum line length 79 characters when possible.
- Functions only used inside the module shouldn't be defined in the
  module header. They can be defined in the beginning of the .c file if
  necessary.
- Symbols defined in headers prefixed with uvg_ or UVG_.
- Includes in alphabetic order.
