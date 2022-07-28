#!/bin/bash

# conda activate mapping

for FILE in *sam
do

if [[ $FILE =~ (.*)\.mapped\.sam$ ]]; then
  SampleID=${BASH_REMATCH[1]}
fi

# clean up read pairing information and flag with SAMtools:
samtools sort -n -O sam ${FILE} | samtools fixmate -m -O bam - ${SampleID}.mapped.fixmate.bam

# sort the bam-file into coordinate order:
samtools sort -O bam -o ${SampleID}.sorted.bam ${SampleID}.mapped.fixmate.bam

# mark duplicated
samtools markdup -S ${SampleID}.sorted.bam ${SampleID}.sorted.duplicates.bam

done

# conda deactivate mapping

