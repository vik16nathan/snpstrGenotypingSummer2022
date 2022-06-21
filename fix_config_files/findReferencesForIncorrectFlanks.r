library(stringr)
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


list_of_samples=c("1232","1393","1791","2006092","2006174","217","2520669","2599208","34896","34966")

#Read in the primer
#Take the primer as an input
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=1) {
  stop("Must supply primer only!", call.=FALSE)
} 

primer <- args[1]

#Load in file containing "incorrect" flanks
incorrect_flank_file <- as.data.frame(read_tsv(paste(
    "Primer",primer,"_possible_incorrect_flanks.tsv",sep="")))

#Load in unreversed incorrect 5' flanks and reversed incorrect 3' flanks
incorrect_forward_flanks <- incorrect_flank_file[
    which(incorrect_flank_file[,"F or R"]=="forward"),"Incorrect Flank"]

incorrect_reverse_flanks <- incorrect_flank_file[
    which(incorrect_flank_file[,"F or R"]=="reverse"),
    "Reversed Incorrect Flank"]

forward_repeat_motif <- incorrect_flank_file[1,"Repeat Motif"]
reverse_repeat_motif <- reverse_flanking_sequence(forward_repeat_motif)


#Iterate through all primer/sample pairs' flank options
#Find the sample/line number that corresponds to the incorrect flank
#Find the .fasta lines for "incorrect" flanks

#Search through every line to see if an incorrect forward or reverse
#flank is found, and then remove the flank from a list of flanks that 
#have not been found yet

#Note that we want only the flanks/.fasta sequences that are ACTUALLY
#incorrect - this means that the 4 b.p. before the 5' flank or
#after the 3' flank belong to the repeat motif

actually_incorrect_forward_flanks <- c()
actually_incorrect_reverse_flanks <- c()

actually_incorrect_forward_sample_list <- c()
actually_incorrect_reverse_sample_list <- c()

actually_incorrect_forward_fasta_references <- c()
actually_incorrect_reverse_fasta_references <- c()

#Create a vector of positions of where the repeat motifs are located
#within the .fasta files - this will allow us to backtrack and find "actual"
#flanking sequences

forward_motif_positions <- c()
reverse_motif_positions <- c()

incorrect_forward_flanks_remaining <- incorrect_forward_flanks
incorrect_reverse_flanks_remaining <- incorrect_reverse_flanks

for(sample in list_of_samples)
{
    if(length(incorrect_forward_flanks_remaining)==0 &&
            length(incorrect_reverse_flanks_remaining)==0) {
                break
    }
    
    flank_file <- as.data.frame(read_tsv(
        paste("./separatedFlanks/P-pyrhulla_",sample,
        ".sorted.duplicates_Primer",primer,"_flanks.tsv",sep="")))
    
    fasta_file <- readLines(
        paste("./FastaInputs/P-pyrhulla_",sample,
        ".sorted.duplicates_Primer",primer,".fasta",sep=""))
    
    actual_forward_flanks <- c()

    for(line in c(1:nrow(flank_file))) {
        if(length(incorrect_forward_flanks_remaining)==0 &&
            length(incorrect_reverse_flanks_remaining)==0) {
                break
            }

        if(flank_file[line,"Motif"]==forward_repeat_motif){
            if(flank_file[line,"5'Flank"] %in% incorrect_forward_flanks_remaining) {
                print(paste("Incorrect forward flank found: ",flank_file[line,"5'Flank"]))
                
                incorrect_fasta_reference <- fasta_file[2*line]
                print(paste("Incorrect fasta reference:", incorrect_fasta_reference))
                flank_pos_within_fasta <- str_locate(incorrect_fasta_reference,flank_file[line,"5'Flank"])[1]
                motif_pos_in_fasta <- flank_pos_within_fasta + 
                        incorrect_flank_file[which(incorrect_flank_file[,"Incorrect Flank"]==flank_file[line,"5'Flank"]),"Motif Position"] -1
                
                num_matching_bp <- 0
                window <- substr(incorrect_fasta_reference, motif_pos_in_fasta-4, motif_pos_in_fasta-1)
                for(j in c(1:4))
                {
                    if(substr(window,j,j)==substr(forward_repeat_motif,j,j)){
                        num_matching_bp <- num_matching_bp + 1
                    }
                }
                if(num_matching_bp >= 3){
                        forward_motif_positions <- c(forward_motif_positions,motif_pos_in_fasta)
                        actually_incorrect_forward_flanks <- c(actually_incorrect_forward_flanks, flank_file[line,"5'Flank"])
                        actually_incorrect_forward_fasta_references <- c(actually_incorrect_forward_fasta_references,
                                                                        incorrect_fasta_reference)
                        
                        actually_incorrect_forward_sample_list <- c(actually_incorrect_forward_sample_list, sample)
                        
                    }
                incorrect_forward_flanks_remaining <- incorrect_forward_flanks_remaining[-which(
                incorrect_forward_flanks_remaining==flank_file[line,"5'Flank"])]

            }
        } else{
            if(flank_file[line,"5'Flank"] %in% incorrect_reverse_flanks_remaining) {
                print(paste("Incorrect reverse flank found: ",flank_file[line,"5'Flank"]))
                incorrect_fasta_reference <- fasta_file[2*line]
                print(paste("Incorrect fasta reference:", incorrect_fasta_reference))
                
                flank_pos_within_fasta <- str_locate(incorrect_fasta_reference,flank_file[line,"5'Flank"])[1]
                print(flank_pos_within_fasta)
                motif_pos_in_fasta <- flank_pos_within_fasta + 
                incorrect_flank_file[which(incorrect_flank_file[,"Reversed Incorrect Flank"]==flank_file[line,"5'Flank"]),"Motif Position"] - 1             
                
                num_matching_bp <- 0
                window <- substr(incorrect_fasta_reference, motif_pos_in_fasta+1, motif_pos_in_fasta+4)
                for(j in c(1:4))
                {
                    if(substr(window,j,j)==substr(reverse_repeat_motif,j,j)){
                        num_matching_bp <- num_matching_bp + 1
                    }
                }
                if(num_matching_bp >= 3){
                        reverse_motif_positions <- c(reverse_motif_positions,motif_pos_in_fasta)
                        actually_incorrect_reverse_flanks <- c(actually_incorrect_reverse_flanks, flank_file[line,"5'Flank"])
                        actually_incorrect_reverse_fasta_references <- c(actually_incorrect_reverse_fasta_references,
                                                                        incorrect_fasta_reference)
                        
                        actually_incorrect_reverse_sample_list <- c(actually_incorrect_reverse_sample_list, sample)
                        
                        
                    }
                incorrect_reverse_flanks_remaining <- incorrect_reverse_flanks_remaining[-which(
                incorrect_reverse_flanks_remaining==flank_file[line,"5'Flank"])]
            
            }
        }
    }

}

#Construct the data frame to be outputted
#Store all incorrect flanks in the orientation in which they were
#found in the .fasta files - do not reverse the 5' flanks, but
#reverse the 3' flanks, since they came from reverse reads

actually_incorrect_reverse_flanks_5p_to_3p <- unlist(lapply(actually_incorrect_reverse_flanks,reverse_flanking_sequence))
all_incorrect_flanks <- c(actually_incorrect_forward_flanks,actually_incorrect_reverse_flanks_5p_to_3p)

all_incorrect_flanks_reversed <- c(actually_incorrect_forward_flanks, actually_incorrect_reverse_flanks)

#Create a vector to represent whether a flank is forward or reverse
full_forward_or_reverse <- unlist(c(rep("forward",length(actually_incorrect_forward_flanks)),
                                    rep("reverse",length(actually_incorrect_reverse_flanks))))

#Create a vector to represent where the repeat motif is located within the fasta input
#for all incorrect_flanks
full_incorrect_flank_motif_positions <- unlist(c(forward_motif_positions,reverse_motif_positions))

#Create a vector of all incorrect reference fastas
full_incorrect_fasta <- unlist(c(actually_incorrect_forward_fasta_references,
                                actually_incorrect_reverse_fasta_references))
#Create output data frame
vector_of_repeat_motifs <- rep(forward_repeat_motif, length(all_incorrect_flanks))
output_df <- data.frame(vector_of_repeat_motifs, all_incorrect_flanks,
        all_incorrect_flanks_reversed, full_incorrect_flank_motif_positions,
        full_forward_or_reverse, full_incorrect_fasta)

colnames(output_df) <- c("Repeat Motif","Incorrect Flank",
"Reversed Incorrect Flank","Motif Position","F or R","Fasta Sequence")

write.table(unique(output_df),
        paste("./Primer",primer,"_actual_incorrect_flanks_with_reference.tsv",sep=""),
            sep="\t",row.names=FALSE,quote=FALSE)
