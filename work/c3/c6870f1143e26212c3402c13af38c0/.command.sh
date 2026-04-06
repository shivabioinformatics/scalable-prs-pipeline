#!/bin/bash -ue
python3 "/Users/labuser3/Documents/BINFX410-Note/homework/Final Project/prs_pipeline/bin/store_results.py" \
    --scores prs_scores.tsv \
    --gwas gwas_summary_stats.tsv \
    --qc-stats qc_stats.tsv \
    --db prs_results.db
