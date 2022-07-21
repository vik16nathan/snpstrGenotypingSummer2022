library(readr)

count_num_ref_and_alt_snps <- function(fasta_lines, snp_position, ref, alt) {
    num_ref_snps <- 0
    num_alt_snps <- 0
    for(i in seq(2,length(fasta_lines),2)) {
        read <- fasta_lines[i]
        snp <- as.character(substr(read, snp_position, snp_position))
        if(snp == ref) {
            num_ref_snps <- num_ref_snps + 1
        } else if(snp == alt) {
            num_alt_snps <- num_alt_snps + 1
        } else {

        }
    }
    print(c(num_ref_snps, num_alt_snps))
}
#Load in intermediate table that needs to be filled in
intermediate_table <- read_tsv("P_pyrhulla-SNPs_before_masking.table_no_STR_alleles_no_masking_intermediate.table")
str_position_table <- read.table("Pyrhulla_pyrhulla_PrimerRef_perf_fixed.bed")

#Iterate through all rows and all samples
sample_names <- c("1232","1393","1791","2006092","2006174","217","2520669","2599208","34896","34966")
for(row in c(1:nrow(intermediate_table))) {
    primer_arg_string <- intermediate_table[row,"CHROM"]
    primer <- as.numeric(substr(primer_arg_string,7,nchar(primer_arg_string)))
    str_start_in_ref <- str_position_table[primer,2]
    str_end_in_ref <- str_position_table[primer,3]

    for(sample in sample_names) {

        snp_zygosity <- intermediate_table[row,paste0(sample,"_SNP_zygosity")]
        str_zygosity <- intermediate_table[row,paste0(sample,"_STR_zygosity")]
        if(snp_zygosity == "he" && str_zygosity == "he") {
            #locate allele 1/allele 2 fasta files
            allele_1_filename <- paste0("./FastaInputs/P-pyrhulla_",sample,"_",primer_arg_string,"_allele1.fasta")
            allele_2_filename <- paste0("./FastaInputs/P-pyrhulla_",sample,"_",primer_arg_string,"_allele2.fasta")
            if(file.exists(allele_1_filename)) {
                allele_1_fasta_lines <- readLines(allele_1_filename)
            } else {
                print(paste("Error! Allele 1 .fasta does not exist for primer",primer,"and sample",sample,"!"))
                next
            }
            if(file.exists(allele_2_filename)) {
                allele_2_fasta_lines <- readLines(allele_2_filename)
            } else {
                print(paste("Error! Allele 2 .fasta does not exist for primer",primer,"and sample",sample,"!"))
                next
            }

            #Find the snp position
            snp_position <- intermediate_table[row, "POS"]
            snp_ref <- as.character(intermediate_table[row, "REF"])
            snp_alt <- as.character(intermediate_table[row, "ALT"])
            print(paste("SNP Ref:", snp_ref))
            print(paste("SNP Alt:", snp_alt))
            if(snp_position > str_end_in_ref) {
                #Find the number of b.p. after the STR where the SNP is located. This will help us find the SNP
                #in the read files, which have more/less STR repeats than in the reference genome
                num_bp_after_str <- snp_position - str_end_in_ref

                #Load in the bed files with start/ending positions of STRs 
                allele_1_bed_filename <- paste0("./bedForMasking/P-pyrhulla_",sample,"_",primer_arg_string,"_allele_1.bed")
                allele_2_bed_filename <- paste0("./bedForMasking/P-pyrhulla_",sample,"_",primer_arg_string,"_allele_2.bed")
                
                if(file.exists(allele_1_bed_filename)) {
                    allele_1_bed_file <- as.data.frame(read_tsv(allele_1_bed_filename,col_names=FALSE))
                    if(any(is.na(allele_1_bed_file))) {
                        print("NAs found in Allele 1 bed file!!")
                        next
                    }
                } else {
                    print("Allele 1 .bed file not found")
                    next
                }
                if(file.exists(allele_2_bed_filename)) {
                    allele_2_bed_file <- as.data.frame(read_tsv(allele_2_bed_filename,col_names=FALSE))
                    if(any(is.na(allele_2_bed_file))) {
                        print("NAs found in Allele 2 bed file!!")
                        next
                    }
                } else {
                    print("Allele 2 .bed file not found")
                    next
                }
                #Iterate through the reads and see which SNP allele is most common in each one
                #Start by iterating through allele 1 .fasta file

                #write a function for this
                allele_1_snp_position <- allele_1_bed_file[3] + num_bp_after_str
                allele_1_ref_alt <- count_num_ref_and_alt_snps(allele_1_fasta_lines, allele_1_snp_position, snp_ref, snp_alt)
                num_allele_1_ref_snps <- allele_1_ref_alt[1]
                num_allele_1_alt_snps <- allele_1_ref_alt[2]
                print(paste("Num allele 1 ref SNPs:", num_allele_1_ref_snps))
                print(paste("Num allele 1 alt SNPs:", num_allele_1_alt_snps))

                allele_2_snp_position <- allele_2_bed_file[3] + num_bp_after_str
                allele_2_ref_alt <- count_num_ref_and_alt_snps(allele_2_fasta_lines, allele_2_snp_position, snp_ref, snp_alt)
                num_allele_2_ref_snps <- allele_2_ref_alt[1]
                num_allele_2_alt_snps <- allele_2_ref_alt[2]
                print(paste("Num allele 2 ref SNPs:", num_allele_2_ref_snps))
                print(paste("Num allele 2 alt SNPs:", num_allele_2_alt_snps))

            } else if(snp_position < str_start_in_ref) {

                allele_1_ref_alt <- count_num_ref_and_alt_snps(allele_1_fasta_lines, snp_position, snp_ref, snp_alt)
                num_allele_1_ref_snps <- allele_1_ref_alt[1]
                num_allele_1_alt_snps <- allele_1_ref_alt[2]
                print(paste("Num allele 1 ref SNPs:", num_allele_1_ref_snps))
                print(paste("Num allele 1 alt SNPs:", num_allele_1_alt_snps))

                allele_2_ref_alt <- count_num_ref_and_alt_snps(allele_2_fasta_lines, snp_position, snp_ref, snp_alt)
                num_allele_2_ref_snps <- allele_2_ref_alt[1]
                num_allele_2_alt_snps <- allele_2_ref_alt[2]
                print(paste("Num allele 2 ref SNPs:", num_allele_2_ref_snps))
                print(paste("Num allele 2 alt SNPs:", num_allele_2_alt_snps))
                
            } else {
                print(paste("Error! SNP found within STR at position", snp_position, "for primer", primer,"and sample", sample,"!"))
            }
        }
    }
}