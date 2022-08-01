#!/bin/bash -i

sampleListFile=$1

./separateFastqByAllele.sh $sampleListFile
#need R script findSTRpositionsFastqEachAllele.r

./cutUpPrimerRefFasta.sh
./mergeSeparatedAlleleReads.sh $sampleListFile
./mapMergedAlleleSeparatedFastq.sh $sampleListFile
./mappingPostMergedSeparated.sh $sampleListFile
./CallVariantsMergedSeparated.sh

#delete all empty files
find ./separatedGVCFs/ -type f -empty -print -delete

#BY HAND - need to make a .contigs list of all Primers with non-empty results in
#separatedGVCFs directory from ./CallVariantsMergedSeparated.sh

#Contigs list should contain one primer per line, entries should be in the form Primer1, etc

#findSNPsWithMergingPart2.sh will consolidate, filter, and convert the SNP results into output tables.
