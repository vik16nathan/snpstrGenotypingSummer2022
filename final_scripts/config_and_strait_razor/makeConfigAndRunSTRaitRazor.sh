#!/bin/bash -i

#First, make .config files without using multiple flanks/sliding window
#Use PERF and have each flank be 10 b.p. long
ref=$1
sampleListFile=$2
cp $ref PrimerRef.fasta
PERF -i PrimerRef.fasta #output filename: PrimerRef_perf.tsv
Rscript PERF_tsv_to_STRaitRazor_config.r 10

#Use more complicated process to fix .config files
./makeConfigFromAllReadsFasta.sh $sampleListFile
#need extractFlankingRegions.r and filterFlankingRegions.r

./fixAllConfig.sh 
#need findIncorrectFlanks.r, findReferencesForIncorrectFlanks.r, findTrueFlanks.r, and fixConfigFile.r

./runAllSTRaitRazor.sh $sampleListFile
#need processMultipleFlanksSTRaitRazorResults.r

Rscript getSTRaitRazorGenotypes.r 0.15 $sampleListFile #15% occurrence threshold for a "significant" allele
