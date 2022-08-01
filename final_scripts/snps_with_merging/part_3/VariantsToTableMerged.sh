#!/bin/bash -i

#EDIT THIS LINE BEFORE RUNNING THIS SCRIPT!!
sampleNames=($(sed -n '62,62p' MergedVariantsBeforeFiltering.vcf | awk '{for (i=10; i<NF-1; i++) printf "\""$i"\"" " "; printf "\""$i"\""}'))
fSampleNames="" 
for sampleName in "${sampleNames[@]}"; do fSampleNames="$fSampleNames -F $sampleName "; done

#extract the primer/sample/allele column names from the .vcf file
gatk VariantsToTable \
    -V MergedVariantsBeforeFiltering.vcf \
    -F CHROM -F POS -F REF -F ALT -F FORMAT $fSampleNames \
    -O MergedVariantsBeforeFiltering.vcf.table
