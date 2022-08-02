library(readr)
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=1) {
  stop("Must supply name of .vcf file with header removed!!", call.=FALSE)
} 


columns_to_remove <- c( "ID", "QUAL", "FILTER", "INFO", "FORMAT")
annotated_vcf <- as.data.frame(read_tsv(args[1]))
annotated_vcf <- annotated_vcf[, 
  !colnames(annotated_vcf) %in% columns_to_remove]

#Extract all non-SNPs
non_snp_vcf <- annotated_vcf[(union(which(nchar(as.character(annotated_vcf[,3])) != 1),
            which(nchar(as.character(annotated_vcf[,4])) != 1))),]

#Filter out all longer variations, leaving only the indels
indel_vcf <- non_snp_vcf[(intersect(which(nchar(as.character(non_snp_vcf[,3])) < 4),
            which(nchar(as.character(non_snp_vcf[,4])) < 4))),]

#Filter out all indels detected within the STR regions
ref_bed <- as.data.frame(read.table("PrimerRef_perf_fixed.bed"))
rows_to_remove <- c()
for(row in c(1:nrow(indel_vcf))) {
    primer_string <- indel_vcf[row, "#CHROM"]
    primer <- as.numeric(substr(primer_string, 7,nchar(primer_string)))
    str_start <- ref_bed[primer,2]
    str_end <- ref_bed[primer,3]
    if(indel_vcf[row,"POS"] <= str_end && indel_vcf[row,"POS"] >= str_start) {
        rows_to_remove <- c(rows_to_remove, row)
    }
}

filtered_indel_vcf <- indel_vcf[which(!(c(1:nrow(indel_vcf)) %in% rows_to_remove)),]
print("Filtered Indel VCF:")
print(filtered_indel_vcf)
#Filter out all indels with too few reads of the alternative allele
rows_with_too_few_reads <- c()
for(row in c(1:nrow(filtered_indel_vcf))) {
    #Find the first occurrence of the alternative allele
    primer_string <- filtered_indel_vcf[row, "#CHROM"]
    primer <- substr(primer_string, 7, nchar(primer_string))
    #print(primer)
    #print(nchar(primer))
    #Extract all columns with sample/allele genotypes corresponding to a particular primer
    primer_columns <- which(substr(colnames(filtered_indel_vcf),1,nchar(primer)) == primer)
    #print(primer_columns)
    indel_found <- FALSE
    for(column in primer_columns) {
        genotype_entry <- filtered_indel_vcf[row,column]
        #A SNP is incorrect if the only occurrence is a 0/1 genotype with < 10 "alternative" reads
        split_entry <- unlist(strsplit(genotype_entry, split=":"))
        ref_alt_zero_one <- split_entry[1]
        #print(ref_alt_zero_one)
        #Case 0: no genotype - next
        if(ref_alt_zero_one == ".") {
            next
        }
        #Case 1: two alternative alleles, meaning that the indel that we identified is a proper indel        
 
        if(substr(ref_alt_zero_one,1,1) != "0" && substr(ref_alt_zero_one,3,3) != "0") {
            indel_found <- TRUE
            break
        }
        #Case 2: one alternative allele - need to see if there is a significant number of reads
        if(substr(ref_alt_zero_one,1,1) == "0" && substr(ref_alt_zero_one,3,3) != "0") {
            ref_alt_read_numbers <- split_entry[2]

            #Extract the number of reads with 0 (reference)
            #or with 1 (alt) SNP
            num_ref_reads <- as.numeric(unlist(strsplit(ref_alt_read_numbers,","))[1])
            num_alt_reads <- as.numeric(unlist(strsplit(ref_alt_read_numbers,","))[2])
            print(num_alt_reads)
            #If the number of alternative read is too small, we haven't actually found an indel
            #THRESHOLD OF 10 READS - change this and see if this makes a difference
            read_threshold <- 10
            if(num_alt_reads > read_threshold) {
                indel_found <- TRUE
                break
            } else {
            }
            
        }
    }
    if(!indel_found) {
        rows_with_too_few_reads <- c(rows_with_too_few_reads, row)
    }
}

read_filtered_indel_vcf <- filtered_indel_vcf[which(!(c(1:nrow(filtered_indel_vcf)) %in% rows_with_too_few_reads)),]
print("Read filtered .vcf:")
print(read_filtered_indel_vcf)
