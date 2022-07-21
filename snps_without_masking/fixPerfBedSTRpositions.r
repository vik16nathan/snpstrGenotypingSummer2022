#Before running this script - run PERF -i Pyrhulla_pyrhulla_PrimerRef.fasta > Pyrhulla_pyrhulla_PrimerRef_perf.tsv
library(readr)
p_pyrhulla_perf_tsv <- as.data.frame(read_tsv("Pyrhulla_pyrhulla_PrimerRef_perf.tsv",
                                              col_names=FALSE))
column_of_fours <- rep(4,30)
p_pyrhulla_bed_output_initial <- p_pyrhulla_perf_tsv[,c(1,2,3,7,8)]
p_pyrhulla_bed_output <- cbind(p_pyrhulla_bed_output_initial[,c(1:3)], column_of_fours)
p_pyrhulla_bed_output <- cbind(p_pyrhulla_bed_output,p_pyrhulla_bed_output_initial[,c(4,5)])
#fix off-by-one error in starting positions
p_pyrhulla_bed_output[,2] <- 1+p_pyrhulla_bed_output[,2]
write.table(p_pyrhulla_bed_output, "./Pyrhulla_pyrhulla_PrimerRef_perf_not_fixed.bed",
            row.names=FALSE,col.names=FALSE,quote=FALSE)

#LOOK THROUGH THIS OUTPUT FILE BY HAND! Copy the results into Pyrhulla_pyrhulla_PrimerRef_perf_fixed.bed
#to see if there are any skipped point mutations/indels in the STR