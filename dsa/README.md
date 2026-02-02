# NanoSeq DSA

## Usage

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

## Development

To generate `compile_commands.json`:

```sh
bear -- ./build.sh
```
