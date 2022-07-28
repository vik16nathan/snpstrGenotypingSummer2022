library(data.table)
library(readr)

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2) {
  stop("Must supply primer and sample!", call.=FALSE)
} 

primer <- args[1]
sample <- args[2]

primer_sample_prefix_string <- paste(sample,".sorted.duplicates_Primer", primer,sep="")

strait_razor_unprocessed_results <- as.data.frame(read_tsv(
                paste("results/",primer_sample_prefix_string,
                    "_STRaitRazor_multiple_flanks.txt",sep=""),
                col_names=FALSE))

#Retain only unique alleles
unique_alleles <- unique(strait_razor_unprocessed_results[,c(1:3)])


#Add fifth and sixth columns, which will later contain the number of forward and reverse occurrences 
#per unique allele
unique_alleles <- cbind(unique_alleles, rep(0,nrow(unique_alleles)))
unique_alleles <- cbind(unique_alleles, rep(0,nrow(unique_alleles)))

#Find the number of forward/reverse occurrences per unique allele 
#(multiple lines could correspond to the same allele
#due to multiple forward/reverse flanking regions in the .config file)
#add up the number of forward and reverse occurrences for each occurrence of a
#unique allele

for(i in c(1:nrow(unique_alleles))) {    
   for(j in c(1:nrow(strait_razor_unprocessed_results))) {
    if(all(unique_alleles[i,c(1:3)] ==
        strait_razor_unprocessed_results[j,c(1:3)])) {
            unique_alleles[i,4] <- unique_alleles[i,4] + strait_razor_unprocessed_results[j,4]
            unique_alleles[i,5] <- unique_alleles[i,5] + strait_razor_unprocessed_results[j,5]
    }
   }
}

#sort in reverse order of total number of occurrences
unique_alleles <- unique_alleles[order(unique_alleles[,4]+unique_alleles[,5], decreasing=TRUE),]

write.table(unique_alleles, paste("results/",primer_sample_prefix_string,
                    "_STRaitRazor_multiple_flanks_unique.txt",sep=""), sep="\t",
                    quote=FALSE,row.names=FALSE, col.names=FALSE)
