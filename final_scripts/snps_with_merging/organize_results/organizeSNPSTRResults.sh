#!/bin/bash -i

#This needs to be run after findSNPsWithMergingPart3.sh 
#Before running this, edit this file with the line number where the table column names are 
#in the .vcf file

#For makeFinalOutputTable.r:
#We also need to store the primer, starting position, ending position, and repeat motif 
#of the STR within the reference genome in .bed file format (without any column names, and separated by a space/tab)
#EDIT THESE POSITIONS BY HAND!!
#... this could also make the original, non-fixed .config files a bit less erroneous, since we have to do this by hand anyways

sampleListFile=$1
tableStart=$(grep -n "#CHROM" MergedVariantsFilteredAnnotated.vcf | cut -d : -f 1)
tail -n +$tableStart MergedVariantsFilteredAnnotated.vcf > MergedVariantsFilteredAnnotatedTableOnly.vcf

./VariantsToTable2Merged.sh

Rscript processFinalVariantTable.r MergedVariantsFilteredAnnotated.table MergedVariantsFilteredAnnotatedTableOnly.vcf

mkdir finalExcelOutputs

Rscript findPerfBedSTRpositions.r
#IMPORTANT - BY HAND - make sure the _fixed.bed file has the proper starting and ending positions of the STRs
#PERF isn't always perfect!
cp PrimerRef_perf_not_fixed.bed PrimerRef_perf_fixed.bed
Rscript makeFinalOutputTable.r $sampleListFile