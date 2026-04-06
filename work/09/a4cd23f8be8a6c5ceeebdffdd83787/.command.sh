#!/bin/bash -ue
echo "Imputing missing genotypes..."
echo "Reference panel: reference_panel.vcf"
echo "Note: In production, this runs Beagle with a 1000 Genomes reference panel"

# demo imputation: replace ./. with 0/0 (simplest possible imputation)
# in production Beagle does this properly using LD patterns
sed 's|\./\.|0/0|g' normalized.vcf > imputed.vcf

# count how many sites were imputed
TOTAL=$(grep -v "^#" normalized.vcf | wc -l)
MISSING=$(grep -v "^#" normalized.vcf | grep -c "\./\." || true)
echo "Imputed ${MISSING} missing calls across ${TOTAL} variant sites"
echo ""
echo "Production workflow:"
echo "  1. java -jar beagle.jar gt=<input> ref=<panel> out=<output>"
echo "  2. Filter: bcftools filter -i 'INFO/DR2 > 0.3' <output>"
