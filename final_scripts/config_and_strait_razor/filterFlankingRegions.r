library(data.table)
library(readr)
#Define a function to convert the 3' to 5' reverse read to a
#5' to 3' forward read (useful for reverse flanking regions)
reverse_flanking_sequence <- function(sequence) {

    #First, reverse the sequence
    splits <- strsplit(sequence, "")[[1]]
    reversed_vector <- rev(splits)

    #Then, convert A <--> T and C <--> G
    converted_bp_vector <- c()
    for(i in c(1:length(reversed_vector))) {
        old_bp <- reversed_vector[i]
        new_bp <- ""
        if(old_bp == "A") {
            new_bp <- "T"
        } else if(old_bp == "T") {
            new_bp <- "A"
        } else if(old_bp == "C") {
            new_bp <- "G"
        } else if(old_bp == "G") {
            new_bp <- "C"
        } else{
            print("Invalid base pair")
        }
        converted_bp_vector <- c(converted_bp_vector, new_bp)
    }
    reversed_sequence <- paste(converted_bp_vector, collapse = "")
    print(reversed_sequence)
}
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2) {
  stop("Must supply primer and list of samples!!", call.=FALSE)
} 


primer <- args[1]
print("Primer:")
print(primer)
print("****************************************************")

list_of_samples <- readLines(args[2])

#Determine forward and reverse directions with reference genome PERF output
perf_reference <- as.data.frame(read_tsv(
            "PrimerRef_perf.tsv",col_names=FALSE))
forward_motif <- perf_reference[primer,8]
    
significant_forward_flanks <- c()
significant_reverse_flanks <- c()

for(sample in list_of_samples)
{
    print("Sample:")
    print(sample)
    primer_sample_prefix_string <- paste(
                sample,".sorted.duplicates_Primer",
                primer,sep="")
    
    flank_filename <- paste("./separatedFlanks/",
                            primer_sample_prefix_string, "_flanks_sorted.tsv",sep="")
    if(!file.exists(flank_filename))
    {
        print("Flanks file does not exist")
        next
    }
    flanks_sorted <- na.omit(as.data.frame(read_tsv(flank_filename,
                            quoted_na=TRUE)))

    flanks_sorted["Frequency"] <- as.numeric(sapply(strsplit(
                            flanks_sorted[,1]," "), `[`, 1))


    flanks_sorted[,1]<- sapply(strsplit(
                            flanks_sorted[,1]," "), `[`, 2)
    #Extract the total number of reads
    total_num_reads <- sum(flanks_sorted["Frequency"])
    if(total_num_reads < 10)
    {
        print(paste("Error: too few reads for primer", primer, 
                "and sample", sample))

    }
    else {
        #Extract the two most frequent flanking region lines in forward
        #and reverse directions (at most four lines)

        #Extract a line only if it has >7.5% of reads
        #(this may seem arbitrary, but we had a 15% threshold for reads
        #for significant alleles after merging forward and reverse reads
        #with STRait razor)

        i <- 1
        forward_flanks <- c()
        while(length(forward_flanks) < 2 && i <= nrow(flanks_sorted))
        {
            if(flanks_sorted[i,"Motif"]==forward_motif &&
                flanks_sorted[i,"Frequency"] >= 0.075*total_num_reads) {
                forward_flanks <- c(forward_flanks, flanks_sorted[i,"5'Flank"])
            }
            i <- (i + 1)
        }

        i <- 1
        reverse_flanks <- c()
        while(length(reverse_flanks) < 2 && i <= nrow(flanks_sorted))
        {
            if(flanks_sorted[i,"Motif"]!=forward_motif &&
                flanks_sorted[i,"Frequency"] >= 0.075*total_num_reads) {
                reversed_flanking_sequence <- reverse_flanking_sequence(
                                                flanks_sorted[i,"5'Flank"])
                reverse_flanks <- c(reverse_flanks, reversed_flanking_sequence)
            }
            i <- (i + 1)
        }

        forward_flanks <- unique(forward_flanks)
        print("Forward flanks:")
        print(forward_flanks)
        reverse_flanks <- unique(reverse_flanks)
        print("Reverse flanks:")
        print(reverse_flanks)
        if(length(forward_flanks)>0){
            significant_forward_flanks <- unique(c(significant_forward_flanks,
                                                     forward_flanks))
        }
        if(length(reverse_flanks)>0){ 
            significant_reverse_flanks <- unique(c(significant_reverse_flanks,
                                                    reverse_flanks))
        }
    }

}

all_flank_pairs <- expand.grid(significant_forward_flanks,
                significant_reverse_flanks)

#If we couldn't find any new flanking regions, use flanks from reference genome
#see PERF_tsv_to_STRaitRazor_config.r
reference_config_df <- as.data.frame(read_tsv(paste("./config/Primer",
                            primer,"_flank10.config",sep="")))
if(nrow(all_flank_pairs) == 0)
{
    output_df <- reference_config_df
} else {
    #Create output config file
    config_file_header <- c("#Marker", "Type", 
            "5'Flank", "3'Flank",	"Motif",	"Period",	"Offset")

    #create columns
    all_markers <- rep(forward_motif, nrow(all_flank_pairs)) #same as motif column
    all_types <- rep("AUTOSOMAL", nrow(all_flank_pairs))
    all_periods <- rep(4, nrow(all_flank_pairs))
    all_offsets <- rep(0, nrow(all_flank_pairs))


    output_df <- data.frame(all_markers, all_types, 
                    all_flank_pairs, all_markers, all_periods, all_offsets)

    colnames(output_df) <- config_file_header
    output_df <- unique(rbind(output_df, reference_config_df))
}
write.table(output_df, 
            paste("config/Primer",primer,
            "_multiple_flanks.config",sep=""), 
                                row.names=FALSE, sep="\t",quote=FALSE)
