library(data.table)

#Start by reading in a perf file generated from the reference genomes of all thirty primers
#This file contains the start and stop positions of STRs within the reference genomes.
#We will use this information in addition to the .fasta file itself to extract flanking sequences, 
#which are required for the .config file of STRaitRazor.

#Perf input generated using the following command
#PERF -i Pyrhulla_pyrhulla_PrimerRef.fasta
#PERF -i Pyrhulla_pyrhulla_PrimerRef_Probe_.fasta
perf_input <- as.data.frame(read.table("Pyrhulla_pyrhulla_PrimerRef_perf.tsv")) #could have off-by-one

fasta_sequences <- readLines("Pyrhulla_pyrhulla_PrimerRef.fasta")

#Create a config file header
config_file_header <- c("#Marker", "Type", "5'Flank", "3'Flank",	"Motif",	"Period",	"Offset")
#Iterate through all primers and find the flanking sequences for the STR
for(primer in 1:30)
{
    #First, scan the reference genome.

    #Isolate the starting and ending points of the STR
    start_point <- perf_input[primer,2]
    end_point <- perf_input[primer,3]

    #Make flanking regions 10 bp long
    #10 bp is enough to decisively locate the STRs; we don't want it any longer because
    #when the number of reads for a certain primer (e.g. 14) is small, if there is more than one mutation in the read's flanking regions,
    #many samples will simply be discarded, preventing genotyping for certain primers
     

    #Every odd-numbered line contains a label for a primer (e.g. Primer1, Primer 3, etc.)
    #Every even-numbered line contains a reference sequence for a primer.

    left_flank <- c(substr(fasta_sequences[2*primer],start_point-10, start_point))
    right_flank <- c(substr(fasta_sequences[2*primer],end_point+1, end_point+11))

    #prepare output .config file for each primer
    motif <- c(noquote(format(perf_input[primer,8])))
    marker <- c(noquote(motif))
    period <- c(4)
    offset <- c(0)
    type <- c(noquote("AUTOSOMAL"))

    output_df <- data.frame(marker, type, left_flank, right_flank, motif, period, offset)
    colnames(output_df) <- config_file_header

    #Scan the other samples' genomes to see if there are indels that we need to include in the .config file
    
    write.table(output_df, paste("./STRaitRazorGenotyping/config/Pyrhulla_pyrhulla_",perf_input[primer,1],".config",sep=""), 
                                row.names=FALSE, sep="\t",quote=FALSE)
                
    
}