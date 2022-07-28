library(data.table)

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2) {
  stop("Must supply primer and sample!", call.=FALSE)
} 
primer <- args[1]
sample <- args[2]

primer_sample_prefix_string <- paste0(sample,".sorted.duplicates_Primer", primer)

perf_input <- as.data.frame(read.table(paste("./PERFallReads/",
            primer_sample_prefix_string,
            "_perf.tsv",sep="")))

fasta_sequences <- readLines(paste("./FastaInputs/",
                            primer_sample_prefix_string,
                            ".fasta",sep=""))

all_left_flanks <- c()
all_right_flanks <- c()
all_motifs <- c()
output_colnames <- c("5'Flank", "3'Flank","Motif")
for(row in 1:nrow(perf_input))
{
    
    #Isolate the starting and ending points of the STR
    start_point <- perf_input[row,2]
    end_point <- perf_input[row,3]
    
    #Make flanking regions 10 bp long
    #10 bp is enough to decisively locate the STRs; we don't want it any longer because
    #when the number of reads for a certain primer (e.g. 14) is small, if there is more than one mutation in the read's flanking regions,
    #many samples will simply be discarded, preventing genotyping for certain primers

    #note: investigate ideal number later
    left_flank <- c(substr(fasta_sequences[2*row],start_point-9, start_point))
    right_flank <- c(substr(fasta_sequences[2*row],end_point+1, end_point+10))

    #prepare output .config file for each primer
    motif <- c(noquote(format(perf_input[row,8])))

    #append to existing lists
    all_left_flanks <- c(all_left_flanks, left_flank)
    all_right_flanks <- c(all_right_flanks, right_flank)
    all_motifs <- c(all_motifs, motif)
}

output_df <- data.frame(all_left_flanks, all_right_flanks, all_motifs)
colnames(output_df) <- output_colnames
write.table(output_df, paste0("./separatedFlanks/",primer_sample_prefix_string,"_flanks.tsv"),
            row.names=FALSE,quote=FALSE,sep="\t")



