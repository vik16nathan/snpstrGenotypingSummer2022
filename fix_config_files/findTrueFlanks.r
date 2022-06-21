library(readr)

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

#Read in the primer
#Take the primer as an input
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=1) {
  stop("Must supply primer only!", call.=FALSE)
} 

primer <- args[1]

#Load in the output from Step 2:
#actually incorrect flanks and corresponding fasta reads

incorrect_flanks_with_ref <- as.data.frame(read_tsv(
    paste("Primer",primer,"_actual_incorrect_flanks_with_reference.tsv",sep="")))

repeat_motif <- incorrect_flanks_with_ref[1,"Repeat Motif"]
reverse_repeat_motif <- reverse_flanking_sequence(repeat_motif)

#Iterate through each incorrect flank reference .fasta read
#Use a 4 bp sliding window to find the correct flank
#For each incorrect flank/reference:
#Locate flank substring within reference - already in output from Part 2
#move the 4 b.p. sliding window by 4 b.p. at a time
#While window matches the repeat motif with >=3 perfectly matching base pairs:
#Slide window back 4 b.p. until true flanking motif is found

true_flanks <- c()

for(row in c(1:nrow(incorrect_flanks_with_ref))) {

    incorrect_flank_pos <- incorrect_flanks_with_ref[row,"Motif Position"]
    fasta_ref <- incorrect_flanks_with_ref[row,"Fasta Sequence"]
    #Slide 4 b.p. window backwards if we're dealing with a 5' flank and 
    #forwards if we're dealing with a 3' flank

    keep_sliding<-TRUE
    if(incorrect_flanks_with_ref[row,"F or R"]=="forward") {

        while(keep_sliding) {
            num_matching_bp <- 0
            window <- substr(fasta_ref, (incorrect_flank_pos-4), (incorrect_flank_pos-1))
            print(paste("Window:",window))
            for(j in c(1:4))
            {
                if(substr(window,j,j)==substr(repeat_motif,j,j)){
                    num_matching_bp <- num_matching_bp + 1
                }
            }
            if(num_matching_bp < 3)
            {
                keep_sliding <- FALSE
            }
            else {
                incorrect_flank_pos <- incorrect_flank_pos-4
            }
        }
        print("True flank:")
        print(substr(fasta_ref, (incorrect_flank_pos-10),(incorrect_flank_pos-1)))
        true_flanks <- c(true_flanks,(substr(fasta_ref, (incorrect_flank_pos-10),(incorrect_flank_pos-1))))

    } else {
         while(keep_sliding) {
            num_matching_bp <- 0
            window <- substr(fasta_ref, (incorrect_flank_pos-4), (incorrect_flank_pos-1))
            for(j in c(1:4))
            {
                if(substr(window,j,j)==substr(reverse_repeat_motif,j,j)){
                    num_matching_bp <- num_matching_bp + 1
                }
            }
            if(num_matching_bp < 3)
            {
                keep_sliding <- FALSE
            }
            else {
                incorrect_flank_pos <- incorrect_flank_pos+4
            }
        }
        print("True flank:")
        print(reverse_flanking_sequence(substr(fasta_ref, (incorrect_flank_pos-10),(incorrect_flank_pos-1))))
        true_flanks <- c(true_flanks,reverse_flanking_sequence(substr(fasta_ref, (incorrect_flank_pos-10),(incorrect_flank_pos-1))))
    }
}

output_df <- data.frame(incorrect_flanks_with_ref[,c(1:2)],true_flanks,incorrect_flanks_with_ref[,"F or R"])
colnames(output_df) <- c("Repeat Motif","Incorrect Flank","Corrected Flank","F or R")
write.table(output_df, paste("Primer",primer,"_corrected_flanks.tsv",sep=""),sep="\t",
            quote=FALSE,row.names=FALSE)





