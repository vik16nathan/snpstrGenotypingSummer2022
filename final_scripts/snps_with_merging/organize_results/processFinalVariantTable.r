library(readr)
library(stringr)

#Function to determine whether a genotype entry with 0/1 in the .vcf file
#is homozygous or heterozygous for a SNP
determine_het_0_1_genotype <- function(genotype_entry) {
  split_entry <- unlist(strsplit(genotype_entry, split=":"))
  ref_alt_read_numbers <- split_entry[2]

  #Extract the number of reads with 0 (reference)
  #or with 1 (alt) SNP
  num_ref_reads <- as.numeric(unlist(strsplit(ref_alt_read_numbers,","))[1])
  num_alt_reads <- as.numeric(unlist(strsplit(ref_alt_read_numbers,","))[2])
  total_num_reads <- as.numeric(split_entry[3])
  print(total_num_reads)
  if(total_num_reads == 0 || split_entry[3] == ".") { return("ho 0") }

  if((num_ref_reads/total_num_reads) >= 0.15 &&
    (num_alt_reads/total_num_reads) >= 0.15) {
      return("he")
    } else if((num_ref_reads/total_num_reads) >= 0.15) {
      return("ho 0")
    } else {
      return("ho 1")
    }
}

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=3) {
  stop("Must supply name of SNP table after filtering and annotation, 
    , the name of the annotated .vcf file, and the list of samples!!", call.=FALSE)
} 


final_variant_table <- as.data.frame(read.table(args[1],header=TRUE))

#Note: the vcf input has to have all the initial header lines removed
#and contain the tabular data only. This can be easily done with sed/cut/tail/etc

annotated_vcf <- as.data.frame(read_tsv(args[2]))
#Remove all columns before primer/sample/allele

columns_to_remove <- c( "ID", "QUAL", "FILTER", "INFO", "FORMAT")
final_variant_table <- final_variant_table[,
        !colnames(final_variant_table) %in% columns_to_remove]

#IMPORTANT - remove quotes from column names so we can locate entries later!
colnames(final_variant_table) <- noquote(colnames(final_variant_table))

annotated_vcf <- annotated_vcf[, 
  !colnames(annotated_vcf) %in% columns_to_remove]

#Keep only first four columns with primer, location, ref, and alt
#Only keeps SNPs and indels (no STRs)

#No inner join needed - the vcf and variant table will always have the same
#number of columns

output_table <- final_variant_table[,c(1:4)]

#Add forty more columns corresponding to each sample.
#For each sample, we will add three columns - 
#1. STR_zygosity (ripped from STRait Razor)
#2. SNP zygosity (to be computed using .vcf)
#3. STR_alleles_where_SNP_is_found (1; 2; 1,2; or 0)
#4. overall_zygosity; "ho" if SNP and STR are "ho"; 
#"he" otherwise

#Make column names
sample_names <- readLines(args[3])
new_column_names <- c()
for(sample_name in sample_names) {
  for(suffix in c("_STR_zygosity","_SNP_zygosity", 
  "_overall_zygosity", "_STR_alleles_where_SNP_is_found")) {
    new_column_names <- c(new_column_names, paste0(sample_name,suffix))
  }
}

#load in STR information
str_table <- as.data.frame(read_tsv("strait_razor_genotypes_0.15_multiple_flanks.tsv"))

snp_allele_information <- as.data.frame(matrix(0,nrow(output_table),length(new_column_names)))
colnames(snp_allele_information) <- new_column_names
output_table <- cbind(output_table, snp_allele_information)



for(primer in c(1:30)) {
  
  primer_string <- paste0("Primer", primer)
  primer_genotype_rows <- which(annotated_vcf["#CHROM"]==primer_string)
  for(sample in sample_names) {
    #See if STR is homozygous, heterozygous, or not in STRait Razor result file
    if(length(which(str_table["Primer"]==as.numeric(primer) & 
          str_table["Sample"]==as.numeric(sample)))==0) {
          print(paste("STR not in result file for primer",primer,"and sample",sample))
          next
    } else {
      print(paste("Primer:",primer))
      print(paste("Sample:",sample))
      print("******************************************************************************")
      str_row <- which((str_table["Primer"])==as.numeric(primer) & 
                      (str_table["Sample"])==as.numeric(sample))
      
      if(str_table[str_row,"Zygosity"] == "ho") {
        #Iterate through all SNP genotypes and see if the sample is homozygous or heterozygous
        #for the particular SNP at the particular primer

        for(row in primer_genotype_rows) {
          output_table[row, paste0(sample,"_STR_zygosity")] <- str_table[str_row,"Zygosity"]
          #Look at allele 1 of the SNP genotype only
          #ISSUE - gets rid of indels (which often have two alternative alleles)
          #work this out later

          if(length(unlist(strsplit(as.character(annotated_vcf[row,"ALT"]),","))) > 1) {
            print(paste("Too many alternative alleles for entry"))
            next
          }
          entry <- format(annotated_vcf[row, which(colnames(annotated_vcf) == paste(primer, sample, "1",sep="_"))])
          entry_zygosity <- determine_het_0_1_genotype(entry)
          if(entry_zygosity == "ho 0") {
            output_table[row, paste0(sample,"_SNP_zygosity")] <- "ho"
            output_table[row, paste0(sample,"_STR_alleles_where_SNP_is_found")] <- "0"
          } else if(entry_zygosity == "ho 1") {
            output_table[row, paste0(sample,"_SNP_zygosity")] <- "ho"
            output_table[row, paste0(sample,"_STR_alleles_where_SNP_is_found")] <- "1"
          } else {
            output_table[row, paste0(sample,"_SNP_zygosity")] <- "he"
            output_table[row, paste0(sample,"_STR_alleles_where_SNP_is_found")] <- "1"
          }
        }

    } else {
        #Iterate through all SNP genotypes and see if the sample is homozygous or heterozygous
        #for the particular SNP at the particular primer
        for(row in primer_genotype_rows) {
          output_table[row, paste0(sample,"_STR_zygosity")] <- str_table[str_row,"Zygosity"]

          #Look at allele 1 and 2 of the SNP genotype
          allele_1_entry <- annotated_vcf[row, which(colnames(annotated_vcf) == paste(primer, sample, "1",sep="_"))]
          allele_1_zygosity <- determine_het_0_1_genotype(allele_1_entry)

          allele_2_entry <- annotated_vcf[row, which(colnames(annotated_vcf) == paste(primer, sample, "2",sep="_"))]
          allele_2_zygosity <- determine_het_0_1_genotype(allele_2_entry)

          if(allele_1_zygosity %in% c("ho 1","he") && allele_2_zygosity %in% c("ho 1","he")) {
            output_table[row, paste0(sample,"_SNP_zygosity")] <- "ho"
            output_table[row, paste0(sample,"_STR_alleles_where_SNP_is_found")] <- "1,2"

          } else if(allele_1_zygosity == "ho 0" && allele_2_zygosity == "ho 0") {
            output_table[row, paste0(sample,"_SNP_zygosity")] <- "ho"
            output_table[row, paste0(sample,"_STR_alleles_where_SNP_is_found")] <- "0"
          } else if(allele_1_zygosity == "ho 0" && allele_2_zygosity %in% c("ho 1","he")) {
            output_table[row, paste0(sample,"_SNP_zygosity")] <- "he"
            output_table[row, paste0(sample,"_STR_alleles_where_SNP_is_found")] <- "2"
          } else {
            output_table[row, paste0(sample,"_SNP_zygosity")] <- "he"
            output_table[row, paste0(sample,"_STR_alleles_where_SNP_is_found")] <- "1"
          }
        }
      }
    }
  }
}

#IMPORTANT - we want SNPs only! This step will get rid of all indels. Perhaps
#consider changing this step if you want to work with indels later on.
output_table <- output_table[(intersect(which(nchar(as.character(output_table[,3])) == 1),
            which(nchar(as.character(output_table[,4])) == 1))),]

#Iterate through the output table and compute the overall zygosity
for(sample in sample_names) {
  for(row in c(1:nrow(output_table))) {
    if(output_table[row, paste0(sample,"_SNP_zygosity")] == "ho" &
       output_table[row, paste0(sample,"_STR_zygosity")] == "ho") {
        output_table[row, paste0(sample,"_overall_zygosity")] <- "ho"
       } else {
          output_table[row, paste0(sample,"_overall_zygosity")] <- "he"
       }
  }
}
head(output_table)
snp_output_table <- output_table[(intersect(which(nchar(as.character(output_table[,3])) == 1),
            which(nchar(as.character(output_table[,4])) == 1))),]

write.table(snp_output_table, paste(args[1],"_no_STR_alleles_intermediate.table",sep=""),sep="\t",
          quote=FALSE,row.names=FALSE)

