#!/bin/bash -i

#We already used GATK to create a .bam file for every pairing of one specific primer and one specific sample
#This was done in createFastqAllPrimersAllSamples.sh, and the results are stored in the FastqInputs directory

#to be run from STRaitRazor genotyping directory
#Now, we want fasta files for each combination of one primer and one sample
#And we then want to combine all files that correspond to the same primer
listOfSamples=("1232" "1393" "1791" "2006092" "2006174" "217" "2520669" "2599208" "34896" "34966")

#Read in the 
for primer in {1..30}
#for primer in {25..26}
do
    for sample in "${listOfSamples[@]}"
    do
        primerString="Primer"
        primerArgString="${primerString}${primer}"
        primerSamplePrefix="P-pyrhulla_${sample}.sorted.duplicates_${primerArgString}"
        #Create a fasta file for each sample/primer pairing
        
        samtools fasta ./separatedFilteredBams/$primerSamplePrefix.bam > \
        ./FastaInputs/$primerSamplePrefix.fasta

        
        #Use PERF on every primer/sample pairing to get the starting and ending position of each STR for every read
        #Make sure that we only do this for non-empty fasta files (aka primer/sample pairings with enough reads)
        if [  -s ./FastaInputs/$primerSamplePrefix.fasta ]; then

            PERF -i ./FastaInputs/${primerSamplePrefix}.fasta -o ./PERFallReads/${primerSamplePrefix}_perf.tsv
        
            #Use R to get the 10 b.p. before and after each STR for every read, still separated for each primer/sample pair
            Rscript extractFlankingRegions.r $primer $sample
            
            #Use uniq -c to find the most frequent flanking regions 
            (head -n 1 ./separatedFlanks/${primerSamplePrefix}_flanks.tsv && tail -n +2 ./separatedFlanks/${primerSamplePrefix}_flanks.tsv | sort | uniq -c | sort -ru) > ./separatedFlanks/${primerSamplePrefix}_flanks_sorted.tsv

    
        fi        
    done
    #Combine results for each primer
    Rscript filterFlankingRegions.r $primer 

        #Place the most signficant flanking regions in the .config file
    
done


