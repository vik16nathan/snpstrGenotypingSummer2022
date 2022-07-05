library(readr)

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2) {
  stop("Must supply name of SNP table after filtering and annotation, 
    and the name of the annotated .vcf file!", call.=FALSE)
} 


final_variant_table <- as.data.frame(read.table(args[1],header=TRUE))
annotated_vcf <- as.data.frame(read_tsv(args[2]))
#Keep only first four columns with primer, location, ref, and alt
#Only keeps SNPs and indels (no STRs)

#!!!!!
#need to perform an inner join to line up rows of vcf and table
#since vcf hasn't removed STRs yet

output_table <- final_variant_table[,c(1:4)]
#Add ten more columns corresponding to each sample.
#For each sample, we will determine whether allele 1, allele 2,
#both, or neither contain the SNP corresponding to the row.

#Make the sample names a command-line argument later.
sample_names <- c("1232","1393","1791","2006092","2006174","217","2520669","2599208","34896","34966")
snp_allele_positions <- as.data.frame(matrix(0,nrow(output_table),length(sample_names)))
colnames(snp_allele_positions) <- sample_names
output_table <- cbind(output_table, snp_allele_positions)

#Separate the primer, sample, and allele names corresponding to each column 
#so that we can locate columns corresponding to primer/sample/allele
primer_sample_allele_df <- as.data.frame(matrix(unlist(strsplit(colnames(annotated_vcf),"_")),ncol=3,byrow=T))
colnames(primer_sample_allele_df) <- c("Primer","Sample","Allele")
#Iterate through the .vcf file to find alleles
for(row in c(1:nrow(annotated_vcf))) {

  primer_number <- substr(annotated_vcf[row, "#CHROM"], 
  7, nchar(annotated_vcf[row,"#CHROM"]))

  columns_with_genotype <- which(primer_sample_allele_df["Primer"]==primer_number)
  #account for the fact that the columns with genotypes in the .vcf file only start at position 10,
  #but primer_sample_allele_df starts containing "actual" primer/sample/allele data at row 3
  #since there are nine columns that don't contain '_'

  #fix this code later. this random 7 is making me uncomfy

  for(column in columns_with_genotype) {

    sample_number <- format(primer_sample_allele_df[column,"Sample"])
    allele_number <- format(primer_sample_allele_df[column,"Allele"])
    if(substr(annotated_vcf[row,7+column],1,3) != "0/0") {
      if(output_table[row,sample_number]==0) {
        output_table[row,sample_number] <- allele_number
      } else {
         output_table[row,sample_number] <- paste(output_table[row,sample_number],
                                          ",",allele_number,sep="")
      }
     
    }
  }
}

output_table <- output_table[(intersect(which(nchar(as.character(output_table[,3])) <= 3),
            which(nchar(as.character(output_table[,4])) <= 3))),]
write.table(output_table, paste(args[1],"_no_STR_alleles.table",sep=""),sep="\t",
          quote=FALSE,row.names=FALSE)