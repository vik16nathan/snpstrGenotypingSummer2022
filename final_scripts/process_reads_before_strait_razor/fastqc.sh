#!/bin/bash

#conda activate qc

mkdir FastQC

for FILE in *R1.trimmed.fastq
do

if [[ $FILE =~ (.*)_(.*)\.trimmed\.fastq$ ]]; then
  SampleID=${BASH_REMATCH[1]}
  Read=${BASH_REMATCH[2]}

## Quality check with FastQC
fastqc -o FastQC ${FILE} ${SampleID}_R2.trimmed.fastq

fi

done

# conda deactivate qc
