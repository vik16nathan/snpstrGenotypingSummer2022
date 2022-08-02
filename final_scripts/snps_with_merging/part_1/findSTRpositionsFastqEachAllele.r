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
    paste("./config/Primer",primer,"_multiple_flanks.config",sep="")))
 
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


#Load in fastq file
fastq_file <- readLines(
    paste("./FastqInputs/",sample,
    ".sorted.duplicates_Primer",primer,".fastq",sep=""))

all_allele_1_R1_lines <- c()
all_allele_1_R2_lines <- c()
#all_allele_1_R2_reversed_lines <- c()


all_allele_2_R1_lines <- c()
all_allele_2_R2_lines <- c()
#all_allele_2_R2_reversed_lines <- c() 

#Important: allele 2 is always longer than allele 1, so always look for
#allele 2 before looking for allele 1, because allele 1 may 
#always be contained within allele 2

#We want new fastq files separated by alleles so that we can have a consistent 
#STR length in each file, allowing CallVariants to detect only the variability in the 
#SNP genotypes.

#Prepare the new output .fastq files.
#Change the header of the separated fastq files to show
#the primer, sample, and allele

allele_1_R1_output_file_name <- paste("./FastqInputs/",
                sample,"_Primer",primer,"_allele1_R1.fastq",sep="")

allele_2_R1_output_file_name <- paste("./FastqInputs/",
                sample,"_Primer",primer,"_allele2_R1.fastq",sep="")

allele_1_R2_output_file_name <- paste("./FastqInputs/",
                sample,"_Primer",primer,"_allele1_R2.fastq",sep="")

allele_2_R2_output_file_name <- paste("./FastqInputs/",
                sample,"_Primer",primer,"_allele2_R2.fastq",sep="")

#Count the number of allele 1 lines, allele 2 lines, and stutter products
allele_1_num <- 0
allele_2_num <- 0
stutter_num <- 0
num_too_small <- 0
for(i in seq(1,(length(fastq_file)-3),4)) {

    #Extract all four lines that comprise the .fastq information for a particular read.
    #1. id_line
    #2. fasta_line
    #3. direction_line (+/-)
    #4. quality_line
    #fasta_line contains the DNA sequence in which we are looking for the STR
    id_line <- fastq_file[i]
    fasta_line <- fastq_file[i+1]
    direction_line <- fastq_file[i+2]
    quality_line <- fastq_file[i+3]

    original_fasta_line <- fasta_line
    reverse_fasta_line <- reverse_sequence(fasta_line)
    #Use the fasta headers to locate whether we have a forward or reverse read.
    forward_read <- FALSE
    if(!is.na(str_locate(fasta_line, allele_1)[2]) || (zygosity == "he" && !is.na(str_locate(fasta_line, allele_2)[2]))) {
        forward_read <- TRUE
    } else if(!is.na(str_locate(reverse_fasta_line, allele_1)[2]) || (zygosity == "he" && !is.na(str_locate(reverse_fasta_line, allele_2)[2]))) {
        forward_read <- FALSE
    } else {
        num_too_small <- num_too_small + 1
        next
    }

    if(!forward_read){
        fasta_line <- reverse_fasta_line
    }
    #account for the possibility that we are dealing with a reverse read by switching the read from
    #3'-5' to 5'-3' and handling it the same way as we would handle a forward read

    #We only need to look for the 5'-3' orientation of the flanks once we've reversed all the reverse 
    #reads

    if(zygosity == "he") {

        if(!is.na(str_locate(fasta_line, allele_2)[2])) {

            
            #Make sure flanking regions are correct and that we aren't dealing with a stutter product
            #Discard any reads that do not have an EXACT match with one of the flanking regions
            #At some point, make this robust to a small amount of mutation

            #Note off-by-one and that we're using 9 and 0 rather than 10 and 1
            #This is because the indices used for masking need to be one b.p.
            #further than the actual locations (since they aren't masked)

            allele_2_start <- str_locate(fasta_line, allele_2)[1] - 1
            allele_2_end <- str_locate(fasta_line, allele_2)[2]
            fasta_forward_flank <- substr(fasta_line, allele_2_start - 9, 
                                        allele_2_start)
            fasta_reverse_flank <- substr(fasta_line, allele_2_end + 1, 
                                        allele_2_end + 10)
            
            
            if(fasta_forward_flank %in% all_forward_flanks && 
                fasta_reverse_flank %in% all_reverse_flanks) {
                        allele_2_num <- allele_2_num + 1
                        
                } else { stutter_num <- stutter_num + 1 }
            if(forward_read) {
                all_allele_2_R1_lines <- c(all_allele_2_R1_lines, id_line, fasta_line, direction_line, quality_line)
            } else {
                all_allele_2_R2_lines <- c(all_allele_2_R2_lines, id_line, original_fasta_line, direction_line, quality_line)

            }
            next
        } 
    }

    if(!is.na(str_locate(fasta_line, allele_1)[2])) {

            allele_1_start <- str_locate(fasta_line, allele_1)[1] - 1
            allele_1_end <- str_locate(fasta_line, allele_1)[2]
            fasta_forward_flank <- substr(fasta_line, allele_1_start - 9, 
                                        allele_1_start)
            fasta_reverse_flank <- substr(fasta_line, allele_1_end + 1, 
                                        allele_1_end + 10)

            if(fasta_forward_flank %in% all_forward_flanks && 
                fasta_reverse_flank %in% all_reverse_flanks) {
                        allele_1_num <- allele_1_num + 1
            }
             if(forward_read) {
                all_allele_1_R1_lines <- c(all_allele_1_R1_lines, id_line, fasta_line, direction_line, quality_line)
            } else {
                
                all_allele_1_R2_lines <- c(all_allele_1_R2_lines, id_line, original_fasta_line, direction_line, quality_line)
                
            }
        } else {
            num_too_small <- num_too_small + 1
        }
}

print(paste("Number of allele 1 lines:", allele_1_num))
print(paste("Number of allele 2 lines:", allele_2_num))
print(paste("Number of stutter/flanking mismatch lines:", stutter_num))
print(paste("Number of lines where neither allele could be located:", num_too_small))

#write separated fasta files for each allele

writeLines(as.character(all_allele_1_R1_lines), allele_1_R1_output_file_name)
writeLines(as.character(all_allele_1_R2_lines), allele_1_R2_output_file_name)

if(zygosity == "he") {
    writeLines(as.character(all_allele_2_R1_lines), allele_2_R1_output_file_name)
    writeLines(as.character(all_allele_2_R2_lines), allele_2_R2_output_file_name)
} 