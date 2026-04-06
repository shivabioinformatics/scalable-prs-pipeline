#!/bin/bash -ue
echo "Normalizing VCF to hg38..."
echo "Chain file: chain_file.chain"

# for the demo, the data is already normalized
# in production, this would run CrossMap + bcftools norm
cp filtered.vcf normalized.vcf

echo "Normalization complete"
echo "Note: In production, this step would run:"
echo "  1. CrossMap.py vcf <chain> <input> <ref.fa> <output>"
echo "  2. bcftools norm -f <ref.fa> -m -both <input> > <output>"
