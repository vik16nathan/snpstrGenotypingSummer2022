#Before running this script - run PERF -i PrimerRef.fasta > PrimerRef_perf.tsv

library(readr)
perf_tsv <- as.data.frame(read_tsv("PrimerRef_perf.tsv",
                                              col_names=FALSE))
column_of_fours <- rep(4,30)
bed_output_initial <- perf_tsv[,c(1,2,3,7,8)]
bed_output <- cbind(bed_output_initial[,c(1:3)], column_of_fours)
bed_output <- cbind(bed_output,bed_output_initial[,c(4,5)])
#fix off-by-one error in starting positions
bed_output[,2] <- 1+bed_output[,2]
write.table(bed_output, "./PrimerRef_perf_not_fixed.bed",
            row.names=FALSE,col.names=FALSE,quote=FALSE)

#LOOK THROUGH THIS OUTPUT FILE BY HAND! Copy the results into PrimerRef_perf_fixed.bed
#to see if there are any skipped point mutations/indels in the STR