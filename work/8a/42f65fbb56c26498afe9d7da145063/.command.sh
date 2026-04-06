#!/bin/bash -ue
python3 "/Users/labuser3/Documents/BINFX410-Note/homework/Final Project/prs_pipeline/bin/visualize_prs.py" \
    --scores prs_scores.tsv \
    --qc-stats qc_stats.tsv \
    --gwas gwas_summary_stats.tsv \
    --outdir .
