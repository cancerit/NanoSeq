# NanoSeq DSA

## Install

The project may be built with CMake, assuming you have the requisite dependencies
```
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build .
./dsa
```

This project explicitly depends on:
- htslib
- OpenSSL
- libdeflate
- PkgConfig (for finding the other dependencies)

A user-specified htslib install may be provided,
otherwise CMake will attempt to find it automatically.

This project requires >=C++23.

See the CMakeLists.txt file for more details.

## Usage

Minimal:

```sh
# Output table: ./output/dsa.bed.gz
./dsa -I ranges.bed -A x.neat.cram -B x.bam -O ./output/
```

Singularity:

```sh
singularity exec \
    -B "/lustre/...:/extra/:ro" \
    -B "/lustre/...:/data/:ro" \
    -B "/lustre/...:/masks/:ro" \
    -B "/lustre/...:/ref/:ro" \
    -B "/lustre/...:/output/" \
    $SINGULARITY_IMAGE /workspace/dsa \
        -I /extra/ranges.bed \
        -A /data/bulk.bam \
        -B /data/duplex.bam \
        -C /masks/snp.bed \
        -R /ref/human/genome.fa \
        -D /masks/noise.bed \
        -d 2 \
        -Q 30 \
        -M 20 \
        -x 12 \
        -O /output/
```

Command line options:

|Option|Type|Default|Description|
|-|-|-|-|
|`A`|file path|-|Bulk (neat) data (BAM/CRAM)|
|`B`|file path|-|Duplex data (BAM/CRAM)|
|`I`|file path|-|Target genomic ranges (BED)|
|`C`|file path|-|SNP mask (BED)|
|`D`|file path|-|Noise mask (BED)|
|`R`|file path|-|Reference genome sequence (FASTA)|
|`O`|file path|2|DSA table output path|
|`Q`|integer|30|Minimum base quality score|
|`M`|integer|0|Minimum read mapping quality|
|`d`|integer|2|Minimum duplex depth|
|`x`|integer|2|DSA table compression level (1-12)|

Migrating from NanoSeq 3:

- the single target genomic range is replaced by a BED file (`I` option)
- no standard output mode

## Development

The following options are available when calling `cmake ..` to set up
the build for the project:
- To generate `compile_commands.json` add the `-DCOMPILE_COMMANDS=ON` flag to
your `cmake ..` call.
- To compile with full debug information use `-DCMAKE_BUILD_TYPE=Debug -DENABLE_DEBUG_FLAGS=ON`.
- To compile tests use `-DMAKE_TEST=ON`. After building with `cmake --build <...>`, tests
may be ran with `./test-dsa`.

These options may be combined as needed

To build the Docker and Singularity image:

```sh
docker build -t dsa . && rm -f dsa.sif && singularity pull dsa.sif
```
