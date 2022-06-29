library(readr)
library(stringr)

reverse_masked_sequence <- function(sequence) {

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
        } else if(old_bp == "N") {
            new_bp <- "N" #accounts for masking
        } else {
            print("Invalid base pair")
        }
        converted_bp_vector <- c(converted_bp_vector, new_bp)
    }
    reversed_sequence <- paste(converted_bp_vector, collapse = "")
    print(reversed_sequence)
}

#load in input data
#Take the masked fasta file name prefix as an input
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2) {
  stop("Must supply masked fasta file name prefix and forward file name!", call.=FALSE)
} 


prefix <- args[1] #for naming the output file
#Load in the R1 masked fasta file
r1_masked_fasta <- readLines(args[2])
print(r1_masked_fasta)
r2_masked_fasta_lines <- c()
for(i in c(1:length(r1_masked_fasta))){

    if(i %% 2 == 1) { #all odd lines have headers
        r2_masked_fasta_lines <- c(r2_masked_fasta_lines, r1_masked_fasta[i])
    } else {
        r2_masked_fasta_lines <- c(r2_masked_fasta_lines, reverse_masked_sequence(r1_masked_fasta[i]))
    }
}

writeLines(r2_masked_fasta_lines, paste("./maskedFasta/",prefix,"_masked_R2.fasta",sep=""))