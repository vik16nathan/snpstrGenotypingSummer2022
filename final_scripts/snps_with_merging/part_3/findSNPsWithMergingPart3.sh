#!/bin/bash -i
./VariantsToTableMerged.sh
./filteringMerged.sh
./CalcGenotypePostMerged.sh
./FilterGenotypesMerged.sh
./VariantAnnotatorMerged.sh

#IMPORTANT!! look at organizeSNPSTRResults.sh and edit tableStart by hand
#based on the .vcf file
