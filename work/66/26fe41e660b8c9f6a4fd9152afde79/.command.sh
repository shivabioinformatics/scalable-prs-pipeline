#!/bin/bash -ue
python3 "/Users/labuser3/Documents/BINFX410-Note/homework/Final Project/prs_pipeline/bin/qc_filter.py" \
    --vcf sample_input.vcf \
    --out-vcf filtered.vcf \
    --out-stats qc_stats.tsv \
    --min-call-rate 0.95 \
    --min-maf 0.01 \
    --min-depth 10 \
    --min-sample-call-rate 0.90
