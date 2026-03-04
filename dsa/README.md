# NanoSeq DSA

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

## Scripts

### Run within a NanoSeq workflow

Equivalent to the DSA section of `runNanoSeq.py`, used in workflows. Includes validation logic.

Expects the `DSA_EXE` environmnet variable to be set as the path to the DSA executable.

```sh
usage: run_dsa.py [-h] [-j INDEX] [-k MAX_INDEX] [-t THREADS] -R REF -A NORMAL -B DUPLEX [-C SNP] [-D MASK] [-d D] [-q Q] [--no_test] [--dry] [--out OUT] [-v]

options:
  -h, --help            show this help message and exit
  -j, --index INDEX     index of the LSF job array. One based
  -k, --max_index MAX_INDEX
                        maximum index of the LSF job array
  -t, --threads THREADS
                        number of threads (1)
  -R, --ref REF         referene sequence
  -A, --normal NORMAL   normal BAM / CRAM
  -B, --duplex DUPLEX   duplex (tumour) BAM / CRAM
  -C, --snp SNP         SNP BED (gz) file
  -D, --mask MASK       mask BED (gz) file
  -d D                  minimum duplex depth (2)
  -q Q                  minimum base quality for normal (30)
  --no_test             skip BAM format tests, use with caution
  --dry                 print the commands and exit
  --out OUT             path of the output files and scratch directory (.)
  -v, --version         show program's version number and exit
```

### Validate results

Usage:

```sh
./validate_dsa_dir.sh tmpDir/dsa
```

This wraps around a more fine-grained tool that compares a compressed DSA table with a report to check for discrepancies in size and, optionally, MD5.

```sh
usage: validate_dsa.py [-h] [--checksum] dsa report

Validate DSA table against truncation.

positional arguments:
  dsa         Path of dsa.bed.gz file
  report      Path of report.json file

options:
  -h, --help  show this help message and exit
  --checksum
```

## File formats

### Report

Format: JSON.

Expected size and MD5 of the DSA table, compressed and uncompressed, for validation purposes.

*E.g.*:

```json
{
  "output_file": "dsa.bed.gz",
  "uncompressed": {
    "size_bytes": 2745867444,
    "md5": "c6702b58ea1bfc04f852507e1331a989"
  },
  "compressed": {
    "size_bytes": 400656342,
    "md5": "c740e3a732e02b1827b891f46b9712f9"
  },
  "compression_ratio": 0.1459
}
```

## Development

To generate `compile_commands.json`:

```sh
bear -- ./build.sh
```

To build the Docker and Singularity image:

```sh
docker build -t dsa . && rm -f dsa.sif && singularity pull dsa.sif
```
