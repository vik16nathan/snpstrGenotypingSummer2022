#Goal - to combine findIncorrectFlanks.r and findReferencesForIncorrectFlanks.r 
#to find flanks that are actually incorrect 
#and to attach their reference sequences on which sliding window
#can be performed, enabling us to find the "true" flanking regions.

library(readr)
library(stringr)
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

determine_incorrect_flank_and_next_window_position <- function(repeat_motif, extended_flank) {
    
    #Input:  extended flank (10 b.p. of the forward/reverse read + 3 b.p.
    #to ensure that any partially-included repeat motifs within a flank
    #are also corrected)

    #Output: non-negative number with the location of the first b.p.
    #after the repeat motif, if the repeat motif is found (75% match)


    #Note: this isn't the final function needed to determine whether a flank
    #in the .config file needs to be fixed - see 
    #determine_incorrect_flank_from_fasta

    #Create a boolean to represent whether or not we have found
    #the repeat motif in the flanking region
    found <- FALSE
    i <- 1
    while(!found && (i + 4) <= nchar(extended_flank))
    {
        print(substr(extended_flank,i,i+4))
        if(compare_five_bp_with_motif(repeat_motif, 
            substr(extended_flank,i,i+4) != -1)) {
            found <- TRUE
            return(i+compare_five_bp_with_motif(repeat_motif, 
            substr(extended_flank,i,i+4)))
        }
    }
    #edge case: the last four b.p. of the extended flank contain the 
    #repeat motif
    if((i+3) == nchar(extended_flank)) {

        num_matching_bp <- 0
        window <- substr(extended_flank, i, (i+3))
        for(j in c(1:4))
        {
            if(substr(window,j,j)==substr(repeat_motif,j,j)){
                num_matching_bp <- num_matching_bp + 1
            }
        }
        if(num_matching_bp >= 3)
        {
            return(i+4)
        }
    }
    return(-1)

}


insertion_found <- function(five_bp, repeat_motif, list_of_indices) {

    #Helper function to get the number of matches between various b.p. in a five b.p. window
    #and the repeat motif

    #If this number is greater than three, then we have found the repeat motif with
    #the insertion reperesented by list_of_indices
    num_matches <- 0
    for(i in c(1:4)) {
        if(substr(five_bp,list_of_indices[i],list_of_indices[i]) == substr(repeat_motif,i,i)) {
            num_matches <- num_matches + 1
        }
    }
    if(num_matches == 4) {
        return(TRUE)
    } else {
        return(FALSE)
    }
}

deletion_found <- function(three_bp, repeat_motif, list_of_indices) {

    #Helper function to get the number of matches between a 3 b.p. window
    #and three base pairs in the repeat motif at positions denoted by list_of_indices

    #If this number is greater than three, 
    #then we have found the repeat motif with
    #the deletion represented by list_of_indices

    num_matches <- 0
    for(i in c(1:3)) {
        if(substr(three_bp,i,i) == substr(repeat_motif,list_of_indices[i],list_of_indices[i])) {
            num_matches <- num_matches + 1
        }
    }
    if(num_matches == 3) {
        return(TRUE)
    } else {
        return(FALSE)
    }
}
compare_five_bp_with_motif <- function(five_bp, repeat_motif) {

    #Input: Five base pairs that could contain the repeat motif
    #with insertions or deletions.

    #Output: boolean value indicating whether the five base pairs
    #have been found.

    #Case 1: no indels (simplest case)
    num_matching_bp <- 0
    window <- substr(five_bp, 1, 4)
    for(j in c(1:4))
    {
        if(substr(window,j,j)==substr(repeat_motif,j,j)){
            num_matching_bp <- num_matching_bp + 1
        }
    }
    if(num_matching_bp >= 3) {
        return(4)
    } #Note: this also takes care of case 3d: deletion of 
    #b.p. 4 of motif (since there will be a 3/4 match nonetheless)
    
    #Case 2: insertion
    #Case 2a: insertion before motif
    #Case 2b: insertion after b.p. 1 of motif
    #Case 2c: insertion after b.p. 2 of motif
    #Case 2d: insertion after b.p. 3 of motif
    if(insertion_found(five_bp, repeat_motif, c(2:5)) ||
        insertion_found(five_bp, repeat_motif, c(1,3,4,5)) ||
        insertion_found(five_bp, repeat_motif, c(1,2,4,5)) || 
        insertion_found(five_bp, repeat_motif, c(1,2,3,5))) {
            return(5)
    }

    #Case 3: deletion
    #Case 3a: deletion at b.p. 1 of motif
    #Case 3b: deletion at b.p. 2
    #Case 3c: deletion at b.p. 3
    #Case 3d: already handled
    three_bp <- substr(five_bp,1,3) 
    if( deletion_found(three_bp, repeat_motif, c(2:4)) ||
       deletion_found(three_bp, repeat_motif, c(1,3,4)) ||
       deletion_found(three_bp, repeat_motif, c(1,2,4)) ) {
        return(3)
    }
    return(-1)

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
forward_repeat_motif <- config_file[1,1]
reverse_repeat_motif <- reverse_flanking_sequence(forward_repeat_motif)

#Extract forward/reverse flanks to be evaluated using .fasta lines.
#Make sure to reorient the reverse flanks to be in the 3'-5' direction.

forward_flanks <- config_file[,3]

#note!! reverse flanks are now in 3'-5' direction
reverse_flanks <- lapply(config_file[,4],reverse_flanking_sequence)

#Extract "extended" forward/reverse flanks (three extra base pairs)
#to see if our flanking region has a partially-embedded STR in it

list_of_samples <- c("1232","1393","1791","2006092","2006174",
                    "217","2520669","2599208","34896","34966")

forward_flanks_remaining <- unique(forward_flanks)
reverse_flanks_remaining <- unique(reverse_flanks)


incorrect_forward_flanks <- c()
incorrect_forward_flank_positions <- c()

incorrect_reverse_flanks <- c()
incorrect_reverse_flank_positions <- c()

#To find .fasta files for all flanking regions for a particular primer,
#we must first iterate through all samples.
for(sample in list_of_samples)
{
    #Keep iterating until all forward/reverse flanks have been classified
    #as incorrect or correct
    if(length(forward_flanks_remaining)==0 &&
            length(reverse_flanks_remaining)==0) {
                break
    }
    
    #extract the list of flanks/fasta files to find
    #a .fasta reference corresponding to a particular flank
    flank_file <- as.data.frame(read_tsv(
        paste("./separatedFlanks/P-pyrhulla_",sample,
        ".sorted.duplicates_Primer",primer,"_flanks.tsv",sep="")))
    
    fasta_file <- readLines(
        paste("./FastaInputs/P-pyrhulla_",sample,
        ".sorted.duplicates_Primer",primer,".fasta",sep=""))
    
    
    for(line in c(1:nrow(flank_file))) {

        if(flank_file[line,"Motif"]==forward_repeat_motif) {
            if(flank_file[line,"5'Flank"] %in% forward_flanks_remaining) {
                print(paste("Forward flank found: ",flank_file[line,"5'Flank"])) 
                incorrect_fasta_reference <- fasta_file[2*line]
                print(paste("Forward fasta reference:", incorrect_fasta_reference))
                flank_pos_within_fasta <- str_locate(incorrect_fasta_reference,flank_file[line,"5'Flank"])[1]
                
                #take three b.p. after the flank to account for incomplete repeat motif inclusions 
                extended_forward_flank <- substr(incorrect_fasta_reference, 
                                flank_pos_within_fasta, flank_pos_within_fasta+12)
                next_window_pos <- determine_incorrect_flank_and_next_window_position(
                    forward_repeat_motif, extended_forward_flank)
                print(paste("Next window pos:", next_window_pos))
                if(next_window_pos != -1) {
                    
                    next_five_bp <- substr(incorrect_fasta_reference, next_window_pos, next_window_pos+4)
                    if(compare_five_bp_with_motif(next_five_bp, forward_repeat_motif)) {
                        print("Truly incorrect forward flank found after looking through fasta:")
                        print(flank_file[line,"5'Flank"])
                        
                        incorrect_forward_flanks <- c(incorrect_forward_flanks, extended_forward_flank)
                        incorrect_forward_flank_positions <- c(incorrect_forward_flank_positions, 
                            flank_pos_within_fasta+next_window_pos-1)
                    }

                    #Note that we could have skipped a previous deletion, since a deletion of b.p. #4
                    #in the motif is viewed as a 3/4 match with a 4 b.p. region
                    next_five_bp_minus_one <- substr(incorrect_fasta_reference, next_window_pos-1, next_window_pos+3) 
                     if(compare_five_bp_with_motif(next_five_bp_minus_one, forward_repeat_motif)) {
                        print("Truly incorrect forward flank found after looking through fasta:")
                        print(flank_file[line,"5'Flank"])
                        
                        incorrect_forward_flanks <- c(incorrect_forward_flanks, extended_forward_flank)
                        incorrect_forward_flank_positions <- c(incorrect_forward_flank_positions, flank_pos_within_fasta+next_window_pos-2)
                    }
                }
                 #Remove forward flank once we've determined whether it's incorrect or not
                 forward_flanks_remaining <- forward_flanks_remaining[-which(
                                    forward_flanks_remaining==flank_file[line,"5'Flank"])]
            }

           

        } else { #repeat what we did above for reverse flanks
            if(flank_file[line,"5'Flank"] %in% reverse_flanks_remaining) {
                print(paste("Reverse flank found: ",flank_file[line,"5'Flank"]))
                incorrect_fasta_reference <- fasta_file[2*line]
                print(paste("Reverse fasta reference:", incorrect_fasta_reference))
                flank_pos_within_fasta <- str_locate(incorrect_fasta_reference,flank_file[line,"5'Flank"])[1]
                
                #take three b.p. after the flank to account for incomplete repeat motif inclusions 
                extended_reverse_flank <- substr(incorrect_fasta_reference, 
                                flank_pos_within_fasta, flank_pos_within_fasta+12)
                next_window_pos <- determine_incorrect_flank_and_next_window_position(
                    reverse_repeat_motif, extended_reverse_flank)
                print(paste("Next window pos:", next_window_pos))
                if(next_window_pos != -1) {
                    next_five_bp <- substr(incorrect_fasta_reference, next_window_pos, next_window_pos+4)
                    if(compare_five_bp_with_motif(next_five_bp, reverse_repeat_motif)) {
                        print("Truly incorrect reverse flank found after looking through fasta:")
                        print(flank_file[line,"5'Flank"])
                        reverse_flanks_remaining <- reverse_flanks_remaining[-which(
                                    reverse_flanks_remaining==flank_file[line,"5'Flank"])]
                        incorrect_reverse_flanks <- c(incorrect_reverse_flanks, extended_reverse_flank)
                        incorrect_reverse_flank_positions <- c(incorrect_reverse_flank_positions, next_window_pos)
                    }
                    #Note that we could have skipped a previous deletion, since a deletion of b.p. #4
                    #in the motif is viewed as a 3/4 match with a 4 b.p. region
                    next_five_bp_minus_one <- substr(incorrect_fasta_reference, next_window_pos-1, next_window_pos+3) 
                     if(compare_five_bp_with_motif(next_five_bp_minus_one, reverse_repeat_motif)) {
                        print("Truly incorrect reverse flank found after looking through fasta:")
                        print(flank_file[line,"5'Flank"])
                        
                        incorrect_forward_flanks <- c(incorrect_reverse_flanks, extended_reverse_flank)
                        incorrect_forward_flank_positions <- c(incorrect_reverse_flank_positions, flank_pos_within_fasta+next_window_pos-2)
                    }
                }
                #Remove reverse flank once we've determined whether it's incorrect or not
                reverse_flanks_remaining <- reverse_flanks_remaining[-which(
                                        reverse_flanks_remaining==flank_file[line,"5'Flank"])]
            }
        }
    }
}

#reverse all 3' flanks, since they came from reverse reads
reversed_reverse_flanks <- lapply(incorrect_reverse_flanks,reverse_flanking_sequence)

#Store the actual sequences (not T/F) of all incorrect flanks
full_incorrect_flanks <- c()
full_incorrect_flanks <- unlist(c(full_incorrect_flanks, 
    incorrect_forward_flanks, reversed_reverse_flanks))

#Store all incorrect flanks in the orientation in which they were
#found in the .fasta files - do not reverse the 5' flanks, but
#reverse the 3' flanks, since they came from reverse reads
full_reversed_incorrect_flanks <- unlist(c(incorrect_forward_flanks, incorrect_reverse_flanks))

#Create a vector to represent whether a flank is forward or reverse
full_forward_or_reverse <- unlist(c(rep("forward",length(incorrect_forward_flank_positions)),
                                    rep("reverse",length(incorrect_reverse_flank_positions))))

full_incorrect_flank_positions <- c(incorrect_forward_flank_positions,
                                    incorrect_reverse_flank_positions)

#Create output data frame
vector_of_repeat_motifs <- rep(forward_repeat_motif, length(full_incorrect_flanks))
output_df <- data.frame(vector_of_repeat_motifs, full_incorrect_flanks,
        full_reversed_incorrect_flanks, full_incorrect_flank_positions,
        full_forward_or_reverse)

colnames(output_df) <- c("Repeat Motif","Incorrect Extended Flank", 
    "Reversed Incorrect Extended Flank","Next Window Position","F or R")


write.table(unique(output_df),
        paste("./actualIncorrectFlanksWithRef/Primer",primer,"_actual_extended_incorrect_flanks_with_reference.tsv",sep=""),
            sep="\t",row.names=FALSE,quote=FALSE)
