library(readr)
library(stringr)

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

strait_razor_genotypes <- as.data.frame(read_tsv(
    "strait_razor_genotypes_0.15_multiple_flanks.tsv"))

#allele 1 is the shorter of the two alleles
#allele 2 is the longer of the two alleles
#if homozygous, then allele 1 == allele 2
allele_1 <- strait_razor_genotypes[
    which(strait_razor_genotypes[,"Primer"]==primer &
    strait_razor_genotypes[,"Sample"]==sample),"Allele 1"]

print(allele_1)

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


#Prepare the new output .fasta files.
#Change the header of the new fasta files to show
#the primer, sample, and allele
new_allele_1_fasta_header <- paste(">Primer",primer,"_sample",
            sample,"_allele_1",sep="")

new_allele_2_fasta_header <- paste(">Primer",primer,"_sample",
            sample,"_allele_2",sep="")

allele_1_output_file_name <- paste("./FastaInputs/P-pyrhulla_",
                sample,".sorted.duplicates_Primer",primer,"_allele1.fasta",sep="")

allele_2_output_file_name <- paste("./FastaInputs/P-pyrhulla_",
                sample,".sorted.duplicates_Primer",primer,"_allele2.fasta",sep="")

for(i in seq(2,length(fasta_file),2)) {

    fasta_line <- fasta_file[i]
    if(zygosity == "he") {
        if(!is.na(str_locate(fasta_line, allele_2)[2])) {
            allele_2_starting_positions <- c(allele_2_starting_positions, 
                    str_locate(fasta_line, allele_2)[1])
            allele_2_ending_positions <- c(allele_2_ending_positions, 
                    str_locate(fasta_line, allele_2)[2])

            all_allele_2_lines <- c(all_allele_2_lines, new_allele_2_fasta_header)
            all_allele_2_lines <- c(all_allele_2_lines, fasta_line)
            next
        } 
    }

    if(!is.na(str_locate(fasta_line, allele_1)[2])) {
            allele_1_starting_positions <- c(allele_1_starting_positions, 
                    str_locate(fasta_line, allele_1)[1])
            allele_1_ending_positions <- c(allele_1_ending_positions, 
                    str_locate(fasta_line, allele_1)[2])
            all_allele_1_lines <- c(all_allele_1_lines, new_allele_1_fasta_header)
            all_allele_1_lines <- c(all_allele_1_lines, fasta_line)
    
        } 


}

print(allele_1_starting_positions)
print(allele_1_ending_positions)
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

print(top_starting_positions)
print(top_ending_positions)
#write output data frame

output_df <- unique(data.frame(rep(repeat_motif,2),top_starting_positions,top_ending_positions))
colnames(output_df) <- c("Repeat Motif","Starting Position","Ending Position")
write.table(output_df, paste("./strPositions/P-pyrhulla_",sample,"_Primer",primer,"_str_positions.tsv",sep=""),
            sep="\t",quote=FALSE,row.names=FALSE)

writeLines(all_allele_1_lines, allele_1_output_file_name)
if(zygosity == "he"){
    writeLines(all_allele_2_lines, allele_2_output_file_name)
}