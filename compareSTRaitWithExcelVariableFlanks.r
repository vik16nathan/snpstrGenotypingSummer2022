library(tidyft)
library(readr)
library(stringr)

for(flank_length in c(10))
{
    print(paste("Flank length:",flank_length))
    for(threshold in c(0.15))
    {


        #load in data, relabel columns to correspond to 
        #excel ("E") or STRait razor ("S")
        excel_genotypes <- as.data.frame(read_tsv("excel_pyrhulla_genotypes.tsv"))
        colnames(excel_genotypes) <- c("Primer", 
                    "Sample", "Allele 1 E", "Allele 2 E", "Zygosity E")

        #Print statements regarding sample sizes
        print("Number of primer/sample configurations in excel table:")
        print(nrow(excel_genotypes))

        strait_razor_genotypes <- as.data.frame(read_tsv(paste("strait_razor_genotypes_",
                                                        threshold,"_flank",flank_length,".tsv",sep="")))
        colnames(strait_razor_genotypes) <- c("Primer", 
                    "Sample", "Allele 1 S", "Allele 2 S", "Zygosity S")
        print("Number of primer/sample configurations in STRait Razor table:")
        print(nrow(strait_razor_genotypes))

        merged_genotype_df <- inner_join(excel_genotypes, 
                strait_razor_genotypes, by=c("Sample", "Primer"))
        print("Number of primer/sample configurations in both:")
        print(nrow(merged_genotype_df))

        print("Number of excel table entries for which STRait razor did not work")
        print(nrow(excel_genotypes)-nrow(merged_genotype_df))

        #Create a new data frame that is the same as the previous merged df
        #with the exception that the shorter allele from excel/STRait Razor is 
        #always the first allele, and the longer allele is the second allele


        for(i in c(1:nrow(merged_genotype_df)))
        {
            #Swap excel alleles, if necessary
            #1. Shorter alleles always come before longer alleles
            #2. If two alleles are the same size, sort them alphabetically as "1" and "2"
            if(nchar(merged_genotype_df[i,"Allele 1 E"]) > nchar(merged_genotype_df[i,"Allele 2 E"]))
            {
                old_value <- merged_genotype_df[i,"Allele 1 E"]
                merged_genotype_df[i,"Allele 1 E"] <- merged_genotype_df[i,"Allele 2 E"]
                merged_genotype_df[i,"Allele 2 E"] <- old_value
            } else if(nchar(merged_genotype_df[i,"Allele 1 E"]) == nchar(merged_genotype_df[i,"Allele 2 E"]) &
                !identical(str_sort(merged_genotype_df[i,c("Allele 1 E","Allele 2 E")]),merged_genotype_df[i,c("Allele 1 E","Allele 2 E")]))
            {
                old_value <- merged_genotype_df[i,"Allele 1 E"]
                merged_genotype_df[i,"Allele 1 E"] <- merged_genotype_df[i,"Allele 2 E"]
                merged_genotype_df[i,"Allele 2 E"] <- old_value
            }else
            {}
            #Swap STRait alleles, if necessary
            if(nchar(merged_genotype_df[i,"Allele 1 S"]) > nchar(merged_genotype_df[i,"Allele 2 S"]))
            {
                old_value <- merged_genotype_df[i,"Allele 1 S"]
                merged_genotype_df[i,"Allele 1 S"] <- merged_genotype_df[i,"Allele 2 S"]
                merged_genotype_df[i,"Allele 2 S"] <- old_value
            }
            else if(nchar(merged_genotype_df[i,"Allele 1 S"]) == nchar(merged_genotype_df[i,"Allele 2 S"]) &
                !identical(str_sort(merged_genotype_df[i,c("Allele 1 S","Allele 2 S")]),merged_genotype_df[i,c("Allele 1 S","Allele 2 S")]))
            {
                old_value <- merged_genotype_df[i,"Allele 1 S"]
                merged_genotype_df[i,"Allele 1 S"] <- merged_genotype_df[i,"Allele 2 S"]
                merged_genotype_df[i,"Allele 2 S"] <- old_value
            }
            
            else {}
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
                setequal(merged_genotype_df[i,c("Allele 1 E","Allele 2 E")],merged_genotype_df[i,c("Allele 1 S","Allele 2 S")]))
            
            #Extract the particular primer/sample combinations for which there 
            #are discrepancies - see how many one-allele mismatches vs. two
            #allele mismatches there are
            one_equal_allele_vector <- c(one_equal_allele_vector,
                sum(merged_genotype_df[i,c("Allele 1 E","Allele 2 E")]==merged_genotype_df[i,c("Allele 1 S","Allele 2 S")]) == 1)
            
            no_equal_allele_vector <- c(no_equal_allele_vector,
                sum(merged_genotype_df[i,c("Allele 1 E","Allele 2 E")]==merged_genotype_df[i,c("Allele 1 S","Allele 2 S")]) == 0)
        }

        #Calculate the proportion of fully equal primer/sample pairs
        print(sum(equal_genotype_vector)/nrow(merged_genotype_df))

        #Calculate the proportion of half-equal primer/sample pairs 
        print(sum(one_equal_allele_vector)/nrow(merged_genotype_df))

        #Calculate the percentage of alleles overall that are equal
        print((sum(merged_genotype_df[,"Allele 1 E"]==merged_genotype_df[,"Allele 1 S"])+ 
                sum(merged_genotype_df[,"Allele 2 E"]==merged_genotype_df[,"Allele 2 S"])))/(2*nrow(merged_genotype_df))


        #Investigate the alleles for which there is one match
        head(merged_genotype_df[which(one_equal_allele_vector),])

        #Investigate the alleles for which there are no matches
        head(merged_genotype_df[which(no_equal_allele_vector),])
        dim(merged_genotype_df[which(no_equal_allele_vector),])

        #Export results for inspection
        write.table(merged_genotype_df[which(no_equal_allele_vector),],
            paste("alleles/full_mismatch_STRait_razor_alleles_",threshold,"_flank",flank_length,".tsv",sep=""),
            row.names=FALSE,quote=FALSE,sep="\t")

        write.table(merged_genotype_df[which(one_equal_allele_vector),],
            paste("alleles/single_mismatch_STRait_razor_alleles_",threshold,"_flank",flank_length,".tsv",sep=""),
            row.names=FALSE,quote=FALSE,sep="\t")

        write.table(merged_genotype_df[c(which(no_equal_allele_vector),
                                        which(one_equal_allele_vector)),],
                paste("alleles/not_double_match_STRait_razor_alleles_",threshold,"_flank",flank_length,".tsv",sep=""),
                    row.names=FALSE,quote=FALSE,sep="\t")
        
        write.table(merged_genotype_df,
                paste("alleles/full_merged_STRait_excel_alleles_",threshold,"_flank",flank_length,".tsv",sep=""),
                    row.names=FALSE,quote=FALSE,sep="\t")
        

    }
}
