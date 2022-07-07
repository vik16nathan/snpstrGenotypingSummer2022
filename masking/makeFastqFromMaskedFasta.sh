#!/bin/bash -i

#Inputs: masked .fasta reads for each allele in R1 and R2 directions
#.fastq header names for masked .fasta file 
#trimmed .fastq file with quality scores for original reads

#Output: a .fastq file with the following format:
#Line 1: @Primer$primer 
#Line 2: masked .fasta read 
#Line 3: + or - (from .fastq file)
#Line 4: quality for each read (from .fastq file)

listOfSamples=("1232" "1393" "1791" "2006092" "2006174" "217" "2520669" "2599208" "34896" "34966")
for sample in "${listOfSamples[@]}"
do
    for direction in R1 R2
    do
        originalFastqFile=./trimmedR1R2fastq/P-pyrhulla_${sample}_${direction}.trimmed.fastq
        for primer in {1..30}
        do
            for allele in 1 2
        
            do
                echo "Sample: $sample"
                echo "Primer: $primer"
                echo "Direction: $direction"
                echo "Allele: $allele"
                fastqHeaderLines=./R1R2headerLists/${sample}_${primer}_${allele}_${direction}.lst
                maskedFasta=./maskedFasta/P-pyrhulla_${sample}_Primer${primer}_allele${allele}_masked_combined_${direction}.fasta
                
                #Iterate through all lines of each .fasta/.fastq file
                numFastaLines=$(wc -l $maskedFasta | cut -d " " -f1)
                if [ $numFastaLines -gt 0 ]
                then
                    (for fastaLineNumber in $(seq 2 2 $numFastaLines)
                    do
                        fastqHeaderLineNumber=$(( fastaLineNumber / 2 ))
                        fastqHeader=$(sed -n "$fastqHeaderLineNumber""p" $fastqHeaderLines)
                        echo "@Primer$primer"
                        #extract the fasta read
                        echo $(sed -n "${fastaLineNumber}p" $maskedFasta)
                        fastqFileLineNumber=$(grep -n "$fastqHeader" $originalFastqFile | cut -d : -f1)
                        #extract the +/- and quality score from fastq file
                        echo $(sed -n "$(( fastqFileLineNumber + 2 ))p" $originalFastqFile)
                        echo $(sed -n "$(( fastqFileLineNumber + 3 ))p" $originalFastqFile)
                    done) > ./maskedFastq/P-pyrhulla_${sample}_Primer${primer}_allele${allele}_masked_combined_${direction}.fastq
                fi
                echo "*******************************************************************************************************"
            done
        done
    done
done
