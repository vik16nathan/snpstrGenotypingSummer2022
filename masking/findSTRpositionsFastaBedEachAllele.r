library(readr)
library(stringr)

reverse_sequence <- function(sequence) {

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
    return(reversed_sequence)
}

#Load in primer and sample
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2) {
  stop("Must supply primer and sample!", call.=FALSE)
} 

primer <- args[1]
sample <- args[2]

#Load in genotype file for each primer/sample pair
config_file <- as.data.frame(read_tsv(
    paste("./config/Pyrhulla_pyrhulla_Primer",primer,"_multiple_flanks.config",sep="")))

repeat_motif <- config_file[1,1]
all_forward_flanks <- unique(config_file[,"5'Flank"])
print("All forward flanks:")
print(all_forward_flanks)

all_reverse_flanks <- unique(config_file[,"3'Flank"])
print("All reverse flanks:")
print(all_reverse_flanks)

strait_razor_genotypes <- as.data.frame(read_tsv(
    "strait_razor_genotypes_0.15_multiple_flanks.tsv"))

#allele 1 is the shorter of the two alleles
#allele 2 is the longer of the two alleles
#if homozygous, then allele 1 == allele 2
allele_1 <- strait_razor_genotypes[
    which(strait_razor_genotypes[,"Primer"]==primer &
    strait_razor_genotypes[,"Sample"]==sample),"Allele 1"]

zygosity <- strait_razor_genotypes[
    which(strait_razor_genotypes[,"Primer"]==primer &
    strait_razor_genotypes[,"Sample"]==sample),"Zygosity"]

if(zygosity == "he") {
    allele_2 <- strait_razor_genotypes[
        which(strait_razor_genotypes[,"Primer"]==primer &
        strait_razor_genotypes[,"Sample"]==sample),"Allele 2"]
}


#Load in fasta file
fasta_file <- readLines(
    paste("./FastaInputs/P-pyrhulla_",sample,
    ".sorted.duplicates_Primer",primer,".fasta",sep=""))

#Iterate through every line of the fasta file
#only even-numbered lines have reads; the odd-numbered lines are headers

#Store all starting/ending positions of STRs
all_allele_1_lines <- c() #store lines for the newly-separated .fasta files by allele
all_allele_2_lines <- c() 
allele_1_starting_positions <- c()
allele_1_ending_positions <- c()

allele_2_starting_positions <- c()
allele_2_ending_positions <- c()

#Important: allele 2 is always longer than allele 1, so always look for
#allele 2 before looking for allele 1, because allele 1 may 
#always be contained within allele 2


#We want new fasta files separated by alleles so that we can mask
#the reads using the starting/ending position of the STR
#and then identify SNPs.

#Prepare the new output .fasta files.
#Change the header of the separated fasta files to show
#the primer, sample, and allele

new_allele_1_fasta_header <- paste("Primer",primer,sep="")

new_allele_2_fasta_header <- paste("Primer",primer, sep="")

allele_1_output_file_name <- paste("./FastaInputs/P-pyrhulla_",
                sample,"_Primer",primer,"_allele1.fasta",sep="")

allele_2_output_file_name <- paste("./FastaInputs/P-pyrhulla_",
                sample,"_Primer",primer,"_allele2.fasta",sep="")

#Count the number of allele 1 lines, allele 2 lines, and stutter products
allele_1_num <- 0
allele_2_num <- 0
stutter_num <- 0
num_too_small <- 0
for(i in seq(2,length(fasta_file),2)) {

    fasta_line <- fasta_file[i]
    reverse_fasta_line <- reverse_sequence(fasta_line)
    #account for the possibility that we are dealing with a reverse read by switching the read from
    #3'-5' to 5'-3' and handling it the same way as we would handle a forward read

    #We only need to look for the 5'-3' orientation of the flanks once we've reversed all the reverse 
    #reads
    if(zygosity == "he") {
        if(!is.na(str_locate(fasta_line, allele_2)[2]) || !is.na(str_locate(reverse_fasta_line, allele_2)[2])) {
            #Make sure flanking regions are correct and that we aren't dealing with a stutter product
            #Discard any reads that do not have an EXACT match with one of the flanking regions
            #At some point, make this robust to a small amount of mutation

            #Note off-by-one and that we're using 11 and 2 rather than 10 and 1
            #This is because the indices used for masking need to be one b.p.
            #further than the actual locations (since they aren't masked)
            if(!is.na(str_locate(reverse_fasta_line, allele_2)[2])){
                fasta_line <- reverse_fasta_line
            }

            allele_2_start <- str_locate(fasta_line, allele_2)[1] - 1
            allele_2_end <- str_locate(fasta_line, allele_2)[2]
            fasta_forward_flank <- substr(fasta_line, allele_2_start - 9, 
                                        allele_2_start)
            fasta_reverse_flank <- substr(fasta_line, allele_2_end + 1, 
                                        allele_2_end + 10)

            if(fasta_forward_flank %in% all_forward_flanks && 
                fasta_reverse_flank %in% all_reverse_flanks) {
                        allele_2_num <- allele_2_num + 1
                        allele_2_starting_positions <- c(allele_2_starting_positions, 
                                                        allele_2_start)
                        allele_2_ending_positions <- c(allele_2_ending_positions, 
                                                        allele_2_end)
                } else { stutter_num <- stutter_num + 1 }

            all_allele_2_lines <- c(all_allele_2_lines, paste(">",new_allele_2_fasta_header,sep=""))
            all_allele_2_lines <- c(all_allele_2_lines, fasta_line)
            next
        } 
    }

    if(!is.na(str_locate(fasta_line, allele_1)[2]) || !is.na(str_locate(reverse_fasta_line, allele_1)[2])) {

            if(!is.na(str_locate(reverse_fasta_line, allele_1)[2])){
                fasta_line <- reverse_fasta_line
            }
            allele_1_start <- str_locate(fasta_line, allele_1)[1] - 1
            allele_1_end <- str_locate(fasta_line, allele_1)[2]
            fasta_forward_flank <- substr(fasta_line, allele_1_start - 9, 
                                        allele_1_start)
            fasta_reverse_flank <- substr(fasta_line, allele_1_end + 1, 
                                        allele_1_end + 10)

            if(fasta_forward_flank %in% all_forward_flanks && 
                fasta_reverse_flank %in% all_reverse_flanks) {
                        allele_1_num <- allele_1_num + 1
                        allele_1_starting_positions <- c(allele_1_starting_positions, 
                                                        allele_1_start)
                        allele_1_ending_positions <- c(allele_1_ending_positions, 
                                                        allele_1_end)
            }

            all_allele_1_lines <- c(all_allele_1_lines, paste(">",new_allele_1_fasta_header,sep=""))
            all_allele_1_lines <- c(all_allele_1_lines, fasta_line)
    
        } else {
            num_too_small <- num_too_small + 1
        }
}

print(paste("Number of allele 1 lines:", allele_1_num))
print(paste("Number of allele 2 lines:", allele_2_num))
print(paste("Number of stutter/flanking mismatch lines:", stutter_num))
print(paste("Number of lines where neither allele could be located:", num_too_small))

#Sort positions from most to least frequent
allele_1_start_frequencies <- as.matrix(sort(table(allele_1_starting_positions), decreasing=TRUE))
allele_1_end_frequencies <- as.matrix(sort(table(allele_1_ending_positions), decreasing=TRUE))

#The row names are the positions; the values are the frequencies
#Choose the first row name in the list (the position with the highest frequency)

top_starting_positions <- c(rownames(allele_1_start_frequencies)[1])
top_ending_positions <- c(rownames(allele_1_end_frequencies)[1])

if(zygosity == "he") {
    allele_2_start_frequencies <- as.matrix(sort(table(allele_2_starting_positions), decreasing=TRUE))
    allele_2_end_frequencies <- as.matrix(sort(table(allele_2_ending_positions), decreasing=TRUE))

    top_starting_positions <- c(top_starting_positions, c(rownames(allele_2_start_frequencies)[1]))
    top_ending_positions <- c(top_ending_positions, c(rownames(allele_2_end_frequencies)[1]))
}


#Consider the first two starting/ending position options
#the rownames are the positions, and the elements are the 
#frequencies for starting/ending positions with frequency vectors
#write separated fasta files for each allele

writeLines(all_allele_1_lines, allele_1_output_file_name)
if(zygosity == "he") {
    writeLines(all_allele_2_lines, allele_2_output_file_name)
}

#Write separated .bed files with starting/ending STR positions for masking later on
if(zygosity == "ho") {
    output_df <- data.frame(new_allele_1_fasta_header,top_starting_positions,top_ending_positions)
    write.table(output_df, paste("./bedForMasking/P-pyrhulla_",sample,"_Primer",primer,"_allele_1.bed",sep=""),
            sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
}
if(zygosity == "he") {
    allele_1_output_df <- data.frame(new_allele_1_fasta_header,top_starting_positions[1],top_ending_positions[1])
    write.table(allele_1_output_df, paste("./bedForMasking/P-pyrhulla_",sample,"_Primer",primer,"_allele_1.bed",sep=""),
            sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
    
    allele_2_output_df <- data.frame(new_allele_2_fasta_header,top_starting_positions[2],top_ending_positions[2])
    write.table(allele_2_output_df, paste("./bedForMasking/P-pyrhulla_",sample,"_Primer",primer,"_allele_2.bed",sep=""),
            sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
}