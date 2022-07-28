#!/bin/bash

#Take in the species name preciding the .fastq extension as an argument
speciesPrefix=$1
for FILE in ${speciesPrefix}-*_R1_*.fastq
do

echo -e "Element:\t${FILE}"

if [[ $FILE =~ (.*-.*)-(.*)_(.*)_(.*)_(.*)_001\.fastq$ ]]; then
  SampleID=${BASH_REMATCH[2]}
  S=${BASH_REMATCH[3]}
  L=${BASH_REMATCH[4]}
  Run=${BASH_REMATCH[5]}

## Trimming with fastp
fastp --correction --cut_right --thread 2 --html ${SampleID}_fastp.html --json ${SampleID}_fastp.json -l 100 -i ${FILE} -I ${speciesPrefix}-${SampleID}_${S}_${L}_R2_001.fastq -o ${SampleID}_R1.trimmed.fastq -O ${SampleID}_R2.trimmed.fastq
# minimum length requirement is specified with -l --length_required --> 100bp

fi
done

echo "done!"
