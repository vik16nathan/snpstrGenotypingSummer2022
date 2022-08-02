#!/bin/bash -i

sampleListFile=$1
contigFile=$2
(./makeConsolidateGVCFs.sh $contigFile) > ConsolidateGVCFsMergedSeparated.sh
chmod +x ConsolidateGVCFsMergedSeparated.sh
./ConsolidateGVCFsMergedSeparated.sh

./JointCallCohortMerged.sh

#IMPORTANT - by hand - determine the line number of the start of the .vcf table (column names) in the gvcf file!!
#choose the number by hand and edit VariantsToTableMerged.sh
gunzip MergedVariantsBeforeFiltering.vcf.gz

#The number needed for VariantsToTableMerged.sh: $(grep -n "#CHROM" MergedVariantsBeforeFiltering.vcf | cut -d : -f 1)
