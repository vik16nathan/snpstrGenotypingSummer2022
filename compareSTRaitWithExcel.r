library(tidyft)
library(readr)

excel_genotypes <- as.data.frame(read_tsv("excel_pyrhulla_genotypes.tsv"))

strait_razor_genotypes <- as.data.frame(read_tsv("strait_razor_genotypes.tsv"))

merged_genotype_df <- inner_join(excel_genotypes, 
        strait_razor_genotypes, by=c("Sample", "Primer"))

colnames(merged_genotype_df) <- c("Sample","Primer","Allele 1 Excel",
    "Allele 2 Excel","Allele 1 STRait","Allele 2 STRait")

#Create a new data frame that is the same as the previous merged df
#with the exception that the shorter allele from excel/STRait Razor is 
#always the first allele, and the longer allele is the second allele


for(i in c(1:nrow(merged_genotype_df)))
{
    #Swap excel alleles, if necessary
    #make sure that all alleles are sorted by size/in alphabetical order
    if(nchar(merged_genotype_df[i,3]) > nchar(merged_genotype_df[i,4]) || 
        sort(merged_genotype_df[i,3:4])!=merged_genotype_df[i,3:4])
    {
        old_value <- merged_genotype_df[i,3]
        merged_genotype_df[i,3] <- merged_genotype_df[i,4]
        merged_genotype_df[i,4] <- old_value
    }
    #Swap STRait alleles, if necessary
    if(nchar(merged_genotype_df[i,5]) > nchar(merged_genotype_df[i,6]) ||
        sort(merged_genotype_df[i,5:6])!=merged_genotype_df[i,5:6])
    {
        old_value <- merged_genotype_df[i,5]
        merged_genotype_df[i,5] <- merged_genotype_df[i,6]
        merged_genotype_df[i,6] <- old_value
    }
}

#Evaluate the percentage of sample/primer pairs 
#for which STRait Razor and the Excel File have the same results

#Create a vector for which the nth entry corresponds to the nth row 
#of the data frame, with boolean values representing whether the 
#Excel and STRait Razor alleles are exactly the same for a particular
#row of the merged data frame.

equal_genotype_vector <- c()
one_equal_allele_vector <- c()
no_equal_allele_vector <- c()
for(i in c(1:nrow(merged_genotype_df)))
{
    equal_genotype_vector <- c(equal_genotype_vector, 
        setequal(merged_genotype_df[i,3:4],merged_genotype_df[i,5:6]))
    
    #Extract the particular primer/sample combinations for which there 
    #are discrepancies - see how many one-allele mismatches vs. two
    #allele mismatches there are
    one_equal_allele_vector <- c(one_equal_allele_vector,
        sum(merged_genotype_df[i,3:4]==merged_genotype_df[i,5:6]) == 1)
    
    no_equal_allele_vector <- c(no_equal_allele_vector,
        sum(merged_genotype_df[i,3:4]==merged_genotype_df[i,5:6]) == 0)
}

#Calculate the proportion of fully equal primer/sample pairs
print(sum(equal_genotype_vector)/nrow(merged_genotype_df))

#Calculate the proportion of half-equal primer/sample pairs 
print(sum(one_equal_allele_vector)/nrow(merged_genotype_df))

#Calculate the percentage of alleles overall that are equal
print((sum(merged_genotype_df[,3]==merged_genotype_df[,5])+ sum(merged_genotype_df[,4]==merged_genotype_df[,6])))/(2*nrow(merged_genotype_df))


#Investigate the alleles for which there is one match
head(merged_genotype_df[which(one_equal_allele_vector),])

#Investigate the alleles for which there are no matches
head(merged_genotype_df[which(no_equal_allele_vector),])
dim(merged_genotype_df[which(no_equal_allele_vector),])

#Export results for inspection
write.xlsx(merged_genotype_df[which(no_equal_allele_vector),],
    "full_mismatch_STRait_razor_alleles.xlsx",
    row.names=FALSE,quote=FALSE)

write.xlsx(merged_genotype_df[which(one_equal_allele_vector),],
    "single_mismatch_STRait_razor_alleles.xlsx",
    row.names=FALSE,quote=FALSE)