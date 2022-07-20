library(readr)
library(stringr)
library(dplyr)

find_mode <- function(x) {
  u <- unique(x)
  tab <- tabulate(match(x, u))
  u[tab == max(tab)]
}

#Load in the intermediate pre-processed output table from processFinalVariantTable.r
intermediate_table <- as.data.frame(read_tsv("P_pyrhulla-SNPs_masked.table_no_STR_alleles_intermediate.table", show_col_types=FALSE))

list_of_samples <- c("1232","1393","1791","2006092","2006174",
                    "217","2520669","2599208","34896","34966")

#Load in the STR genotype table
str_table <- as.data.frame(read_tsv("strait_razor_genotypes_0.15_multiple_flanks.tsv", show_col_types=FALSE))

#Iterate through every single primer
for(primer in c(1:30)) {
#for(primer in c(12)) {
    #Load in the repeat motif
    config_filename <- paste(
    "./config/Pyrhulla_pyrhulla_Primer",primer,
    "_multiple_flanks.config",sep="")

    if(file.exists(config_filename)) {
        config_file <- as.data.frame(read_tsv(config_filename, show_col_types=FALSE))
    } else {
        print(paste("Config file not found for primer", primer))
        next
    }
    

    #Extract the repeat motif
    repeat_motif <- config_file[1,1]

    #See if there are any SNPs corresponding to primer
    primer_string <- paste0("Primer",primer)
    possible_snp_rows <- which(intermediate_table[,"CHROM"]==primer_string)
    #note - "SNP" is used loosely here. Some of the variants that we want to 
    #process are actually indels, but we are calling everything SNPs because we 
    #ultimately want SNPSTRs

    if(length(possible_snp_rows)==0){
        print(format(paste("No SNPs found for Primer", primer)))
    }

    #Prepare two-row pieces of final output table for each primer
    #Each row corresponds to an STR allele, primer, and sample

    primer_sample_str_lengths <- c()
    primer_sample_str_genotypes <- c()
    sample_names <- c()
    repeat_motifs <- rep(repeat_motif, 2)
    
    overall_primer_output_data_frame <- c()
    #Store all STR starting positions - note that this is sometimes incorrectly computed,
    #so we need to use the starting position that's most common among all samples
    #for a particular primer.
    all_str_starting_positions <- c()

    #Iterate through each sample
    for(sample in list_of_samples) {

        #print(paste("Sample:", sample))
        #Create a vector corresponding to the starting and ending positions of the longer STR
        #This will allow us to filter out any "SNPs" within the STR

        longer_STR_allele_positions <- c()

        sample_names <- rep(sample, 2)
        str_row <- which(str_table[,"Primer"]==primer & str_table[,"Sample"]==sample)
        if(length(str_row)==0) {
            print(format(paste("No STR genotyped for Primer", primer, " and sample", sample)))
            next
        }

        primer_sample_str_genotypes <- c(str_table[str_row,"Allele 1"], str_table[str_row, "Allele 2"])

        #Determine STR zygosity
        sample_STR_zygosity <- str_table[str_row,"Zygosity"]

        #Get the bed file for allele 1
        allele_1_bed_file_name <- paste0("./bedForMasking/P-pyrhulla_",sample,"_",primer_string,"_allele_1.bed")
      
        if(!file.exists(allele_1_bed_file_name)) {
            print(paste("Allele 1 .bed file does not exist for primer",primer,"and sample",sample,"!"))
            
        } else {
            primer_sample_allele_1_bed_file <- read_tsv(allele_1_bed_file_name,
                                                        col_names=FALSE, show_col_types = FALSE)
            if(any(is.na(primer_sample_allele_1_bed_file))) {
                print(paste("NAs in allele 1 .bed file for primer", primer, "and sample", sample))
                next
            }  
            str_starting_position <- primer_sample_allele_1_bed_file[2]
            all_str_starting_positions <- c(all_str_starting_positions, str_starting_position)
            longer_STR_allele_positions <- c(primer_sample_allele_1_bed_file[2], primer_sample_allele_1_bed_file[3])
        }

        primer_sample_allele_1_length <- nchar(str_table[str_row, "Allele 1"])
        primer_sample_allele_1_decimal_length <- paste0(floor(primer_sample_allele_1_length/4), ".", 
                                                (primer_sample_allele_1_length %% 4) )
        if(sample_STR_zygosity == "ho") {
            primer_sample_str_lengths <- rep(primer_sample_allele_1_decimal_length, 2)
            
        
        } else {
            #Get the bed file for allele 2
            allele_2_bed_file_name <- paste0("./bedForMasking/P-pyrhulla_",sample,"_",primer_string,"_allele_2.bed")
            if(!file.exists(allele_2_bed_file_name)) {
                print(paste("Allele 2 .bed file does not exist for primer",primer,"and sample", sample,"!"))
            } else {
                     primer_sample_allele_2_bed_file <- read_tsv(allele_2_bed_file_name,
                                                        col_names=FALSE, show_col_types = FALSE)
                if(any(is.na(primer_sample_allele_2_bed_file))) {
                    print(paste("NAs in allele 2 .bed file for primer", primer, "and sample", sample))
                } 
                longer_STR_allele_positions <- c(primer_sample_allele_2_bed_file[2], primer_sample_allele_2_bed_file[3])
            }
       
            primer_sample_allele_2_length <- nchar(str_table[str_row,"Allele 2"])
            primer_sample_allele_2_decimal_length <- paste0(floor(primer_sample_allele_2_length/4), ".", 
                                                    (primer_sample_allele_2_length %% 4) )
            primer_sample_str_lengths <- c(primer_sample_allele_1_decimal_length, primer_sample_allele_2_decimal_length)
        }

        two_row_output_df <- cbind(sample_names, primer_sample_str_genotypes, repeat_motifs, primer_sample_str_lengths)
        colnames(two_row_output_df)[1] <- "Sample Name"
        colnames(two_row_output_df)[2] <- "STR Genotypes"
        colnames(two_row_output_df)[3] <- "Repeat Motif"
        colnames(two_row_output_df)[4] <- "Msat"

        #Iterate through putative SNPs
        overall_zygosity <- "ho"
        if(length(possible_snp_rows) != 0) {
            for(snp_row in possible_snp_rows) {
                #Determine the position of the SNP
                snp_position <- intermediate_table[snp_row, "POS"]
                #see if potential SNP is embedded within STR
                if(length(longer_STR_allele_positions) == 2 && !any(is.na(longer_STR_allele_positions)) &&
                    (snp_position >= longer_STR_allele_positions[1]) && (snp_position <= longer_STR_allele_positions[2])) {
                    print(paste("SNP/indel at position", snp_position, "is within STR"))
                    next
                } else {
                    #If any of the SNPs are heterozygous, then the overall SNPSTR genotype is heterozygous
                    if(intermediate_table[snp_row, paste0(sample,"_overall_zygosity")] == "he") {
                        overall_zygosity <- "he"
                    }
                    #Extract the reference/alternative alelles
                    ref_allele <- intermediate_table[snp_row, "REF"]
                    alt_allele <- intermediate_table[snp_row, "ALT"]
                    
                    #Create a column that will align the SNP alleles with the STRs where they are found (if any)
                    #in the output table
                    snp_position_column <- c()

                    #Look at which STR alleles correspond with the SNP
                    if(intermediate_table[snp_row, paste0(sample, "_STR_alleles_where_SNP_is_found")] == "0") {
                        snp_position_column <- c("-","-")
                    } else if (intermediate_table[snp_row, paste0(sample, "_STR_alleles_where_SNP_is_found")] == "1") {
                        if(intermediate_table[snp_row, paste0(sample, "_SNP_zygosity")] == "ho") {
                            snp_position_column <- c(alt_allele, alt_allele)
                        } else {
                            snp_position_column <- c(alt_allele, "-")
                        }
                    } else if (intermediate_table[snp_row, paste0(sample, "_STR_alleles_where_SNP_is_found")] == "2") {
                        snp_position_column <- c("-", alt_allele)
                    } else {
                        snp_position_column <- c(alt_allele, alt_allele)
                    }
                    two_row_output_df <- cbind(two_row_output_df, snp_position_column)
                    colnames(two_row_output_df)[ncol(two_row_output_df)] <- paste0(snp_position, " ",ref_allele,"->",alt_allele)
                }
            }
        } else { #the overall zygosity is just the STR zygosity in the absence of SNPs
            overall_zygosity <- sample_STR_zygosity
        }
        #Add a column corresponding to the overall zygosity
        two_row_output_df <- cbind(two_row_output_df, rep(overall_zygosity, 2))
        #print(two_row_output_df)
        colnames(two_row_output_df)[ncol(two_row_output_df)] <- "inheritance"

        #append two-column intermediate data frame to overall data frame
        if(length(overall_primer_output_data_frame) == 0) {
            overall_primer_output_data_frame <- two_row_output_df
        } else {
            overall_primer_output_data_frame <- bind_rows(as.data.frame(overall_primer_output_data_frame), 
                                                          as.data.frame(two_row_output_df))
        }
    }

    if(length(overall_primer_output_data_frame) == 0) {
        print(paste("No results for primer",primer))
        next
    }
    #remove all NAs
    overall_primer_output_data_frame <- overall_primer_output_data_frame[, colSums(is.na(overall_primer_output_data_frame))==0]

    #Find the most frequent STR starting position
    most_frequent_STR_starting_position <- as.numeric(find_mode(all_str_starting_positions)) + 1 #correct off-by-one error

    colnames(overall_primer_output_data_frame)[which(colnames(overall_primer_output_data_frame) == "Msat")] <- 
                            paste(most_frequent_STR_starting_position, "Msat")
    
    write.table(overall_primer_output_data_frame, paste0("./finalExcelOutputs/Primer",primer,"_SNPSTRs.tsv"),row.names=FALSE,
        quote=FALSE, sep="\t")

}