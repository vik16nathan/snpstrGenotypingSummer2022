#!/bin/bash

##Reference
Ref=$1
sampleListFile=$2
bwa index $Ref

#Take in the list of samples as a argument
mapfile -t listOfSamples < $sampleListFile
for SampleID in "${listOfSamples[@]}"
do
    FILE=${SampleID}_merged.fastq
    bwa mem -R "@RG\tID:${SampleID}\tLB:${SampleID}\tPL:ILLUMINA\tPM:HISEQ\tSM:${SampleID}"  $Ref ${FILE}  > ${SampleID}.mapped.sam
done