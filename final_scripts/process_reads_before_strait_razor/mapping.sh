#!/bin/bash

##Reference
Ref=$1

bwa index $Ref

for FILE in *R1.trimmed.fastq
do

if [[ $FILE =~ (.*)_(.*)\.trimmed\.fastq$ ]]; then
  SampleID=${BASH_REMATCH[1]}
  Read=${BASH_REMATCH[2]}
fi

bwa mem -R "@RG\tID:${SampleID}\tLB:${SampleID}\tPL:ILLUMINA\tPM:HISEQ\tSM:${SampleID}"  $Ref ${FILE} ${SampleID}_R2.trimmed.fastq  > ${SampleID}.mapped.sam
done
