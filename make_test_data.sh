#!/usr/bin/env bash
# Build a small test dataset: reads overlapping the first 100 kb of the reference + that 100 kb reference.
# NOT yet tested. Needs the full raw data once (data/raw) and the p18 environment.
set -euo pipefail
REF=data/ref/ref.fasta; R=NC_000913.3:1-100000
mkdir -p tests && samtools faidx $REF && samtools faidx $REF $R > tests/ref_100kb.fasta
minimap2 -t 2 -ax map-ont $REF data/raw/ERR14686234.fastq.gz | samtools sort -o /tmp/l.bam - && samtools index /tmp/l.bam
samtools view -b /tmp/l.bam $R | samtools fastq - | gzip > tests/long_100kb.fastq.gz
minimap2 -t 2 -ax sr $REF data/raw/SRR1030394_1.fastq.gz data/raw/SRR1030394_2.fastq.gz | samtools sort -o /tmp/s.bam - && samtools index /tmp/s.bam
samtools view -b /tmp/s.bam $R | samtools sort -n - | samtools fastq -1 tests/short_100kb_R1.fastq.gz -2 tests/short_100kb_R2.fastq.gz -0 /dev/null -s /dev/null -n -
ls -l tests
