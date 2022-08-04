#!/bin/bash -i

#We already used GATK to create a .bam file for every pairing of one specific primer and one specific sample
#This was done in createFastqAllPrimersAllSamples.sh, and the results are stored in the FastqInputs directory

#to be run from STRaitRazor genotyping directory
#Now, we want fasta files for each combination of one primer and one sample
#And we then want to combine all files that correspond to the same primer

#Read in the list of samples
sampleListFile=$1
mapfile -t listOfSamples < $sampleListFile
for primer in {1..30}
do
    for sample in "${listOfSamples[@]}"
    do
        primerString="Primer"
        primerArgString="${primerString}${primer}"
        primerSamplePrefix="${sample}.sorted.duplicates_${primerArgString}"
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
    #Rscript filterFlankingRegions.r $primer $sampleListFile
    Rscript filterFlankingRegionsEvenMore.r $primer $sampleListFile
        #Place the most signficant flanking regions in the .config file
    
done


