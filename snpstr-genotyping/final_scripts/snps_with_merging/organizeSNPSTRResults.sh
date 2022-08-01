#!/bin/bash -i

#This needs to be run after findSNPsWithMergingPart3.sh 
#Before running this, edit this file with the line number where the table column names are 
#in the .vcf file

#IMPORTANT!! Replace 74 with the actual line number where table starts
sampleListFile=$1
tableStart=74
tail -n +$tableStart MergedVariantsFilteredAnnotated.vcf > MergedVariantsFilteredAnnotatedTableOnly.vcf

./VariantsToTable2Merged.sh

Rscript processFinalVariantTable.r MergedVariantsFilteredAnnotated.table MergedVariantsFilteredAnnotatedTableOnly.vcf

mkdir finalExcelOutputs
Rscript makeFinalOutputTable.r $sampleListFile