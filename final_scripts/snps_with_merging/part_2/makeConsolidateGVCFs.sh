#!/bin/bash

#Goal - create ConsolidateGVCFsMaskedSeparated.sh using all g.vcf.gz files in the maskedSeparatedGVCFs directory
#./makeConsolidateGVCFs.sh > ConsolidateGVCFsMergedSeparated.sh

#delete all empty files

contigsFilename=$1
genomicsdb-workspace-path=$2
find ./separatedGVCFs/ -type f -empty -print -delete

echo "#!/bin/bash -i"
echo ""
echo "gatk --java-options \"-Xmx4g -Xms4g\" GenomicsDBImport \\"

for file in $(ls ./separatedGVCFs/ | grep g.vcf.gz$ )
do
    echo "  -V "\"\.\/separatedGVCFs\/$file\" \\""
done

echo "--genomicsdb-workspace-path ${genomicsdb-workspace-path} \\"
echo "-L ${contigsFilename}"
