library(readr)

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=1) {
  stop("Must supply name of SNP table after filtering and annotation!", call.=FALSE)
} 

final_variant_table <- as.data.frame(read.table(args[1],header=TRUE))
#Keep only first four columns with primer, location, ref, and alt
#Only keeps SNPs and indels (no STRs)

print(as.character(final_variant_table[,4]))
output_table <- final_variant_table[(intersect(which(nchar(as.character(final_variant_table[,3])) <= 3),
            which(nchar(as.character(final_variant_table[,4])) <= 3))),c(1:4)]

write.table(output_table, paste(args[1],".noSTR.table",sep=""),
            sep="\t",row.names=FALSE,quote=FALSE)

