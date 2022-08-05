#!/bin/bash -i

sampleListFile=$1

./separateFastqByAllele.sh $sampleListFile
#need R script findSTRpositionsFastqEachAllele.r

./cutUpPrimerRefFasta.sh
./mergeSeparatedAlleleReads.sh $sampleListFile
./mapMergedAlleleSeparatedFastq.sh $sampleListFile
./mappingPostMergedSeparated.sh $sampleListFile
./CallVariantsMergedSeparated.sh $sampleListFile

#delete all empty files
find ./separatedGVCFs/ -type f -empty -print -delete

#BY HAND - need to make a contigs .list file of all Primers in the main working directory.
#Include primers with non-empty results from the separatedGVCFs directory

#Contigs list should contain one primer per line, entries should be in the form Primer1, etc

#findSNPsWithMergingPart2.sh will consolidate, filter, and convert the SNP results into output tables.
#Need to make a contigs .list file of all Primers in the main working directory.
#Include primers with non-empty results from the separatedGVCFs directory
for primer in {1..30}
do
    primerSearchString="Primer$primer""_"
    numGVCFs=$((ls ./separatedGVCFs | grep $primerSearchString ) | wc -l)

    if [ $numGVCFs -gt 0 ]
    then
        echo "Primer$primer" >> contigs.list
    fi

done

contigFile=contigs.list
(./makeConsolidateGVCFs.sh $contigFile) > ConsolidateGVCFsMergedSeparated.sh
chmod +x ConsolidateGVCFsMergedSeparated.sh
#Consolidate all ~300 GVCFs separated by primer/sample/allele into one output directory
./ConsolidateGVCFsMergedSeparated.sh

#Combine all results into one table with all variants for all primer/sample/allele combinations
./JointCallCohortMerged.sh

gunzip MergedVariantsBeforeFiltering.vcf.gz

#Extract the starting point of the .vcf table above and use it to make a table of all variants prior to filtering 
tableStart=$(grep -n "#CHROM" MergedVariantsBeforeFiltering.vcf | cut -d : -f 1)
./VariantsToTableMerged.sh $tableStart

#Filtering
./filteringMerged.sh
./CalcGenotypePostMerged.sh
./FilterGenotypesMerged.sh
#Annotate all final variants after filtering into a .vcf file
./VariantAnnotatorMerged.sh
