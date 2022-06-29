#!/bin/bash

#Goal - create ConsolidateGVCFsMaskedSeparated.sh using all g.vcf.gz files in the maskedSeparatedGVCFs directory
#./makeConsolidateGVCFs.sh > ConsolidateGVCFsMaskedSeparated.sh

#delete all empty files

find ./maskedSeparatedGVCFs/ -type f -empty -print -delete

echo "#!/bin/bash -i"
echo ""
echo "gatk --java-options \"-Xmx4g -Xms4g\" GenomicsDBImport \\"

for file in $(ls ./maskedSeparatedGVCFs/ | grep g.vcf.gz$ )
do
    echo "  -V "\"\.\/maskedSeparatedGVCFs\/$file\" \\""
done

echo "--genomicsdb-workspace-path P-pyrhulla-masked \\"
echo "-L P-pyrhulla_contigs.list"