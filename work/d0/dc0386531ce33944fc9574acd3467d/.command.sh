#!/bin/bash -ue
python3 "/Users/labuser3/Documents/BINFX410-Note/homework/Final Project/prs_pipeline/bin/calculate_prs.py" \
    --vcf imputed.vcf \
    --gwas gwas_summary_stats.tsv \
    --output prs_scores.tsv \
    --p-threshold 0.05
