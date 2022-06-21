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

determine_incorrect_flank <- function(repeat_motif, flank){
    
    #Input:  flank (could be any length)
    #Output: boolean value if flank contains a 4 b.p. region that matches
    #75% with the repeat motif of the STR

    #Create a boolean to represent whether or not we have found
    #the repeat motif in the flanking region
    found <- FALSE
    i <- 1
    while(!found && (i + 3) <= nchar(flank))
    {
        #Take a 4 b.p. sliding window
        #Count the number of complete matches between base pairs
        #in the window and the repeat motif
        num_matching_bp <- 0
        window <- substr(flank, i, (i+3))
        for(j in c(1:4))
        {
            if(substr(window,j,j)==substr(repeat_motif,j,j)){
                num_matching_bp <- num_matching_bp + 1
            }
        }
        if(num_matching_bp >= 3)
        {
            found <- TRUE
        }
        i <- i+1
    }
    print(found) #return true if the flank is INCORRECT

}

determine_motif_position <- function(repeat_motif, incorrect_flank){
    
    #Input:  incorrect flank with 75% match with repeat motif 
    #Output: integer index of where the repeat motif 
    #is located within 10 b.p. flank

    found <- FALSE
    index_found <- 1
    while(!found && (index_found + 3) <= nchar(incorrect_flank))
    {
    
        num_matching_bp <- 0
        window <- substr(incorrect_flank, index_found, (index_found+3))
        for(j in c(1:4))
        {
            if(substr(window,j,j)==substr(repeat_motif,j,j)){
                num_matching_bp <- num_matching_bp + 1
            }
        }
        if(num_matching_bp >= 3)
        {
            found <- TRUE
        }
        else
        {
            index_found <- (index_found + 1)
        }
    }
    print(index_found) #return index

}


#Take the primer as an input
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=1) {
  stop("Must supply primer only!", call.=FALSE)
} 

primer <- args[1]

#Load in the .config file
config_file <- as.data.frame(read_tsv(paste(
    "./config/Pyrhulla_pyrhulla_Primer",primer,
    "_multiple_flanks.config",sep="")))

#Extract the repeat motif
repeat_motif <- config_file[1,1]

#Nomenclature - 5' flank = forward flank, 3' flank = reverse flank
#Evaluate all forward flanking regions
incorrect_forward_flanks <- c() #boolean vector of whether flank is incorrect
for(forward_flank in config_file[,3]){
    incorrect_forward_flanks <- c(incorrect_forward_flanks,
        determine_incorrect_flank(repeat_motif, forward_flank))
}


incorrect_forward_flank_positions <- c()
for(flank in config_file[which(incorrect_forward_flanks),3]){
    
    incorrect_forward_flank_positions <- c(incorrect_forward_flank_positions,
                                determine_motif_position(repeat_motif,flank))
}

incorrect_reverse_flanks <- c()
for(reverse_flank in config_file[,4]){
    incorrect_reverse_flanks <- c(incorrect_reverse_flanks,
        determine_incorrect_flank(repeat_motif, reverse_flank))
}

incorrect_reverse_flank_positions <- c()
for(flank in config_file[which(incorrect_reverse_flanks),4]){
    
    incorrect_reverse_flank_positions <- c(incorrect_reverse_flank_positions,
        determine_motif_position(repeat_motif,flank))
}

#reverse all 3' flanks, since they came from reverse reads
reversed_reverse_flanks <- lapply(config_file[ 
    which(incorrect_reverse_flanks), 4],
                                reverse_flanking_sequence)

#Store the actual sequences (not T/F) of all incorrect flanks
full_incorrect_flanks <- c()
full_incorrect_flanks <- c(full_incorrect_flanks, 
        config_file[which(incorrect_forward_flanks),3])

full_incorrect_flanks <- c(full_incorrect_flanks, 
        config_file[which(incorrect_reverse_flanks),4])

#Store all incorrect flanks in the orientation in which they were
#found in the .fasta files - do not reverse the 5' flanks, but
#reverse the 3' flanks, since they came from reverse reads
full_reversed_incorrect_flanks <- unlist(c(
    config_file[which(incorrect_forward_flanks),3],
                                reversed_reverse_flanks))

#Create a vector to represent whether a flank is forward or reverse
full_forward_or_reverse <- unlist(c(rep("forward",length(incorrect_forward_flank_positions)),
                                    rep("reverse",length(incorrect_reverse_flank_positions))))

full_incorrect_flank_positions <- c(incorrect_forward_flank_positions,
                                    incorrect_reverse_flank_positions)

#Create output data frame
vector_of_repeat_motifs <- rep(repeat_motif, length(full_incorrect_flanks))
output_df <- data.frame(vector_of_repeat_motifs, full_incorrect_flanks,
        full_reversed_incorrect_flanks, full_incorrect_flank_positions,
        full_forward_or_reverse)

colnames(output_df) <- c("Repeat Motif","Incorrect Flank",
"Reversed Incorrect Flank","Motif Position","F or R")

write.table(unique(output_df),
        paste("./Primer",primer,"_possible_incorrect_flanks.tsv",sep=""),
            sep="\t",row.names=FALSE,quote=FALSE)
