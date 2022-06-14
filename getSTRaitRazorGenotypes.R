#Read in all STRaitRazor result files
list_of_samples=c("1232","1393","1791","2006092","2006174","217","2520669","2599208","34896","34966")
list_of_primers_to_include=c(3,4,5,8,10,11,12,13,14,15,16,17,18,19,20,22,23,25,26,28,29,30)

#accept the decimal threshold for identifying an allele as a "significant allele" as a user input
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("Must supply a threshold!", call.=FALSE)
} 
threshold <- args[1]
print(threshold)

#Initialize the columns of the output data frame
all_primers <- c() 
all_samples <- c() 
all_allele_1 <- c()
all_allele_2 <- c()
for(primer in list_of_primers_to_include)
{
    for(sample in list_of_samples)
    {
          #Ensure we're not dealing with a primer/sample configuration for which STRait Razor failed
          result_filename <- paste("./results/P-pyrhulla_",
                                    sample,".sorted.duplicates_Primer",
                                    primer, 
                                    "_STRaitRazor_multiple_flanks_unique.txt",sep="")
          if (!file.exists(result_filename)) {
            print(paste("Strait razor result file does not exist for primer", 
                                primer, "and sample", sample))
            next
          } else if (file.size(result_filename) == 0L) {
            print(paste("Strait razor result file is empty for primer", 
                                primer, "and sample", sample))
            next
          } else {
                    strait_razor_output <- as.data.frame(read.table(result_filename))

                    #Count the total number of reads for a particular primer/sample configuration
                    #total_num_reads <- sum(strait_razor_output[,5])+sum(strait_razor_output[,6])
                    repeat_motif <- substr(strait_razor_output[1,1],1,4)
                    #make sure to count reads only for matches that don't correspond to a flanking region insertion
                    non_inserted_alleles <- which(substr(repeat_motif,1,1) == substr(strait_razor_output[,4],1,1))
                    total_num_reads <- sum(strait_razor_output[non_inserted_alleles,5],
                                            strait_razor_output[non_inserted_alleles,6])
                    #IMPORTANT: reject all primer/sample combinations for which there are <10 reads!!
                    if(total_num_reads >=10)
                    {
                        
                        #Create a vector to represent "significant" alleles, with threshold as a user input
                        significant_alleles <- c()
                        all_primers <- c(all_primers, primer)
                        all_samples <- c(all_samples, sample)
                        for(i in c(1:nrow(strait_razor_output)))
                        {
                            genotype <- format(strait_razor_output[i,4])
                            repeat_motif <- substr(strait_razor_output[i,1],1,4)
                            num_reads_per_allele <- strait_razor_output[i,5]+strait_razor_output[i,6]
                            if(num_reads_per_allele/total_num_reads > threshold && 
                                    substr(genotype,1,1) == substr(repeat_motif,1,1))
                            {
                                significant_alleles <- c(significant_alleles, genotype)
                            }
                        }

                        #IMPORTANT - account for insertions in the flanking region (1E error)
                        #Add significant allele only if its first character is the same as the repeat motif

                        #Note: there is still an edge case that isn't addressed by this - what if there's a point
                        #mutation in the first character of the STR? Then, we won't identify it as a significant
                        #allele, when it is very much significant

                        #note - significant alleles are always in order of frequency
                        if(length(significant_alleles)>0)
                        {
                            significant_alleles <- rev(significant_alleles[order(nchar(significant_alleles), 
                                                                        significant_alleles)])
                            all_allele_1 <- c(all_allele_1,significant_alleles[1])
                            #Case 1 - one significant allele - homozygous
                            if(length(significant_alleles) == 1)
                            {
                                all_allele_2 <- c(all_allele_2,significant_alleles[1])
                            }
                            #Case 2 - two significant alleles  - heterozygous
                            #Case 3 - more than two significant alleles - pick the two most frequent, call it a heterozygote
                            #alleles are already in descending order of frequency based on the ordering of output rows from STRait razor
                            #note that this doesn't happen often - only 1 time out of ~200 configurations
                            else
                            {
                                all_allele_2 <- c(all_allele_2,significant_alleles[2])
                            }
                            all_primers <- c(all_primers, primer)
                            all_samples <- c(all_samples, sample)
                        }
                    }
                   
                }
    }
}

#Create a boolean vector where true represents homozygous and false represents heterozygous
all_homo_hetero <- factor(all_allele_1 == all_allele_2, labels=c("he","ho"))
output_df <- data.frame(all_primers, all_samples, all_allele_1, all_allele_2, all_homo_hetero)
colnames(output_df) <- c("Primer", "Sample", "Allele 1", "Allele 2", "Zygosity")
write.table(output_df, paste("strait_razor_genotypes_",threshold,"_multiple_flanks.tsv", sep=""),
            row.names=FALSE, sep="\t", quote=FALSE)