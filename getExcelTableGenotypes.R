library(readxl)
library(qmrparser)

#Function taken from https://www.geeksforgeeks.org/how-to-read-a-xlsx-file-with-multiple-sheets-in-r/    
multiplesheets <- function(fname) {
   
  # getting info about all excel sheets
  sheets <- readxl::excel_sheets(fname)
  tibble <- lapply(sheets, function(x) readxl::read_excel(fname, sheet = x))
  data_frame <- lapply(tibble, as.data.frame)
  
  # assigning names to data frames
  names(data_frame) <- sheets

  # print data frame
  print(data_frame)
}


#Function to convert the bracket/parenthesis notation in mutated alleles into a full genotype
expand_full_genotype <- function(geno_string)
{
    print(geno_string)
    #Remove all spaces
    geno_string <- format(gsub(" ", "", geno_string, fixed = TRUE))
    #First, expand any parenthesis, if any
    #Start by finding the region within parentheses

    if(grepl("(",geno_string,fixed=TRUE))
    {
        starting_parenthesis_index <- unlist(gregexpr("(",geno_string,fixed=TRUE))
        ending_parenthesis_index <- unlist(gregexpr(")",geno_string,fixed=TRUE))
        #Assume that there is only an integer number of repeats after a parenthesis
        num_parenthesis_repeats <- as.numeric(substring(geno_string, ending_parenthesis_index+1, ending_parenthesis_index+1))
        expanded_parenthesis_string <- ""
        for(i in c(1:num_parenthesis_repeats))
        {
            expanded_parenthesis_string <- paste(expanded_parenthesis_string, substring(geno_string, starting_parenthesis_index+1,
                                                    ending_parenthesis_index-1),sep="")
        }
        geno_string <- paste(substring(geno_string, 1, starting_parenthesis_index-1), expanded_parenthesis_string,
                                substring(geno_string, ending_parenthesis_index+2, nchar(geno_string)),sep="")
        
    }

    #Now, convert the brackets.
    starting_bracket_indices <- as.vector(unlist(gregexpr("[",geno_string,fixed=TRUE)))
    ending_bracket_indices <- as.vector(unlist(gregexpr("]",geno_string,fixed=TRUE)))

    #Iterate through all bracket pairings and expand all except for the last one, which could have a decimal suffix
    length_added <- 0 #keep updating the bracket indices as the intermediate string keeps expanding

    #Figure out if we need to iterate to the second-to-last or the last repetition of []
    #We iterate to the second-to-last repetition if the last character of the genotype string is a number
    #Otherwise, we can iterate to the last repetition
    if(isDigit(substring(geno_string,nchar(geno_string),nchar(geno_string))))
    {
        num_iters <- length(starting_bracket_indices)-1
    }
    else
    {
        num_iters <- length(starting_bracket_indices)
    }
    if(num_iters > 0)
    {
        for(i in c(1:num_iters))
        {
            starting_bracket_index <- starting_bracket_indices[i]+length_added
            ending_bracket_index <- ending_bracket_indices[i]+length_added
            unprocessed_motif <- substring(geno_string, starting_bracket_index+1,
                                                    ending_bracket_index-1)
            #remove all special characters from repeat motif
            repeat_motif <- gsub('[^[:alnum:] ]', '', unprocessed_motif)
            
            num_digits <- 1
            
            if(isDigit(substring(geno_string,ending_bracket_index+2,ending_bracket_index+2)))
            {
                num_digits <- 2
                num_bracket_repeats <- as.numeric(substring(geno_string, ending_bracket_index+1, ending_bracket_index+2))
            }
            else
            {
                num_bracket_repeats <- as.numeric(substring(geno_string, ending_bracket_index+1, ending_bracket_index+1))
            }

            expanded_bracket_string <- ""
            for(i in c(1:num_bracket_repeats))
            {
                expanded_bracket_string <- paste(expanded_bracket_string, repeat_motif,sep="")
            }
            geno_string <- paste(substring(geno_string, 1, starting_bracket_index-1), expanded_bracket_string,
                                    substring(geno_string, ending_bracket_index+1+num_digits, nchar(geno_string)),sep="")
            length_added <- length_added + nchar(repeat_motif)*(num_bracket_repeats-1) - 2 - num_digits
            #subtract each time due to the removal of the brackets themselves + the number of repeats
        }
    }
    
    #First, see if the last character of the genotype string is a number.
    #If so, perform an expansion, with decimals as needed.
    #If not, then do nothing.
    #See if the next-to-last character of the genotype string is a decimal, indicating that we must do a partial expansion

    i <- length(starting_bracket_indices)
    starting_bracket_index <- starting_bracket_indices[i]+length_added
    ending_bracket_index <- ending_bracket_indices[i]+length_added
    if(isDigit(substring(geno_string,nchar(geno_string),nchar(geno_string))))
    {
        final_repeat_number <- as.numeric(substring(geno_string, ending_bracket_index+1, nchar(geno_string)))
        num_digits <- nchar(substring(geno_string, ending_bracket_index+1, nchar(geno_string)))
        num_bracket_repeats <- round(final_repeat_number,0)
        num_extra_chars <- 10*(final_repeat_number %% 1)
        expanded_bracket_string <- ""
        print("Geno String:")
        print(geno_string)
        print("Starting Bracket Index:")
        print(starting_bracket_index)

        print("Ending Bracket Index:")
        print(ending_bracket_index)
        #remove all special characters
        repeat_motif <- gsub('[^[:alnum:] ]', '', 
                            substring(geno_string, starting_bracket_index+1, ending_bracket_index-1))
        for(i in c(1:num_bracket_repeats))
        {
            expanded_bracket_string <- paste(expanded_bracket_string, repeat_motif,sep="")
        }

        geno_string <- paste(substring(geno_string, 1, starting_bracket_index-1), expanded_bracket_string,
                                substring(geno_string, ending_bracket_index+1+num_digits, nchar(geno_string)),sep="")
        
        if(num_extra_chars > 0)
        {
            print("Number of characters of repeat motif:")
            print(nchar(repeat_motif))
            for(j in 1:num_extra_chars)
            {
                geno_string <- paste(geno_string, substring(repeat_motif,j,j),sep="")
            }
            
        } 
    }
    else
    {
        #geno_string <- paste(geno_string, substring(geno_string, ending_bracket_index+1, nchar(geno_string)))
    }
    
    print(geno_string)

}

  
# specifying the path name
path <- "./Results-List-Pyrrhula_pyrrhula_FINAL.xlsx"
all_pyrhulla_dfs <- multiplesheets(path)

#Goal: an output data frame with the following columns:
#1. Primer 
all_primers <- c()
#2. Sample 
all_samples <- c()
#3. Allele 1 (all base pairs)
all_allele_1 <- c()
#4. Allele 2 (all base pairs)
all_allele_2 <- c()

#Iterate through all sheets
list_of_primers_to_include=c(3,4,5,8,10,11,12,13,14,15,16,17,18,19,20,22,23,25,26,28,29,30)

#Create a "dictionary" to go from sample codes in the excel table to the easier numeric values
excel_sample_names <- c("ZFMK-DNA-FD19595294","ZFMK-DNA-FD19595310","ZFMK-DNA-FD19595301","ZFMK-TIS-2006092",
                        "ZFMK-TIS-2006174","ZFMK-DNA-FD19595295","ZFMK-TIS-2520669","ZFMK-TIS-2599208",
                        "ZFMK-TIS-34896","ZFMK-TIS-34966")

numeric_sample_names <- c(1232,1393,1791,2006092,2006174,217,2520669,2599208,34896,34966)
sample_dict_df <- as.data.frame(cbind(excel_sample_names,numeric_sample_names))
colnames(sample_dict_df) <- c("Excel","Numeric")

for(primer in list_of_primers_to_include)
{
    print("Primer:")
    print(primer)
    print("**********************************************")
    #Load in the data frame for a particular primer
    primer_df <- as.data.frame(all_pyrhulla_dfs[paste("PyrPyr",primer,sep="")])
    
    #Locate the column where the microsatellite genotypes are located
    msat_column <- which(grepl("Msat",primer_df[1,]))

    #Only update the columns of the dataframe if the genotype is not NA
    for(i in c(1:10)) #these rows contain the genotypes
    {
        #Row 2*i contains allele 1; row 2*i+1 contains allele 2 corresponding to the same sample
        if(!is.na(primer_df[2*i,msat_column]))
        {
                     
            excel_sample <- primer_df[2*i,1]
            numeric_sample <- sample_dict_df[which(sample_dict_df["Excel"]==excel_sample),"Numeric"]
            all_samples <- c(all_samples, format(numeric_sample))
            all_primers <- c(all_primers, primer)
            #add the repeat motif - later, combine this with the forensic number to get the full genotype
            genotype_1 <- format(primer_df[2*i,5]) 
            #this is the same thing as repeat_motif reference din the function above

            if(nchar(genotype_1) == 4) #case 1: allele 1 is an STR with no mutations 
            {
                full_genotype_1 <- ""
                msat_number_1 <- primer_df[2*i,msat_column]
                #If there are no mutations, then the genotype has at most four characters - get rid of any 
                #mistaken letters after the decimal point
                if(isLetter(substring(msat_number_1,nchar(msat_number_1),nchar(msat_number_1))))
                {
                    msat_number_1 <- substring(msat_number_1,1,nchar(msat_number_1)-1)
                }
                
                print(msat_number_1)
                #Edge case - the third character of either allele is - and not .
                #Replace the - with a .
                if(nchar(msat_number_1)>2 & substring(msat_number_1,3,3) == "-")
                {
                    msat_number_1 <- paste(substring(msat_number_1,1,2),".",
                                                        substring(msat_number_1,4,nchar(msat_number_1)),
                                                        sep="")
                }
                num_repeats <- round(as.numeric(msat_number_1),0)
                num_extra_chars <- (as.numeric(msat_number_1) %% 1)*10
                for(j in c(1:num_repeats))
                {
                    full_genotype_1 <- paste(full_genotype_1, genotype_1, sep="")
                }
                if(num_extra_chars > 0)
                {
                    for(j in 1:num_extra_chars)
                    {
                            full_genotype_1 <- paste(full_genotype_1, substring(genotype_1,j,j),sep="")
                    }
                }
                all_allele_1 <- c(all_allele_1, full_genotype_1)

            } #case 2: allele 1 has mutations - see expand_full_genotype function definition
            #to go from brackets to the genotype

            else{
                all_allele_1 <- c(all_allele_1, expand_full_genotype(genotype_1))
            }

            #repeat both cases for allele 2
            genotype_2 <- format(primer_df[2*i+1,5])
            if(nchar(genotype_2) == 4) #case 1: we have an STR with no mutations
            {
                full_genotype_2 <- ""
                msat_number_2 <- primer_df[2*i+1,msat_column]
                #If there are no mutations, then the genotype has at most four characters - get rid of any 
                #mistaken letters after the decimal point
                if(isLetter(substring(msat_number_2,nchar(msat_number_2),nchar(msat_number_2))))
                {
                    msat_number_2 <- substring(msat_number_2,1,nchar(msat_number_2)-1)
                }
                
                print(msat_number_2)
                 #Edge case - the third character of either allele is - and not .
                 #Replace the - with a .
                if(nchar(msat_number_2)>2 & substring(msat_number_2,3,3) == "-")
                {
                    msat_number_2 <- paste(substring(msat_number_2,1,2),".",
                                                        substring(msat_number_2,4,nchar(msat_number_2)),
                                                        sep="")
                }
                num_repeats <- round(as.numeric(msat_number_2),0)
                num_extra_chars <- (as.numeric(msat_number_2) %% 1)*10
                for(j in c(1:num_repeats))
                {
                    full_genotype_2 <- paste(full_genotype_2, genotype_2, sep="")
                }
                if(num_extra_chars > 0)
                {
                    for(j in 1:num_extra_chars)
                    {
                            full_genotype_2 <- paste(full_genotype_2, substring(genotype_2,j,j),sep="")
                    }
                }
                all_allele_2 <- c(all_allele_2, full_genotype_2)

            }
            else{
                all_allele_2 <- c(all_allele_2, expand_full_genotype(genotype_2))
            }
        }
    }
}

all_homo_hetero <- factor(all_allele_1 == all_allele_2, labels=c("he","ho"))
output_df <- data.frame(all_primers, all_samples, all_allele_1, all_allele_2, all_homo_hetero)
colnames(output_df) <- c("Primer", "Sample", "Allele 1", "Allele 2", "Zygosity")
write.table(output_df, "excel_pyrhulla_genotypes.tsv",row.names=FALSE,quote=FALSE,sep="\t")
