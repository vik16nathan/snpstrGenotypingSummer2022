library(readr)
p_pyrhulla_perf_tsv <- as.data.frame(read_tsv("C:/Users/vik16/OneDrive/Documents/FOGS Summer 2022/Pyrhulla_pyrhulla_PrimerRef_perf.tsv",
                                              col_names=FALSE))
p_pyrhulla_perf_tsv
column_of_fours <- rep(4,30)
p_pyrhulla_bed_output_initial <- p_pyrhulla_perf_tsv[,c(1,2,3,7,8)]
p_pyrhulla_bed_output <- cbind(p_pyrhulla_bed_output_initial[,c(1:3)], column_of_fours)
p_pyrhulla_bed_output <- cbind(p_pyrhulla_bed_output,p_pyrhulla_bed_output_initial[,c(4,5)])
#add one to all starting and ending positions of STRs based on an error discovered
#at the end of the HipSTR pipeline

p_pyrhulla_bed_output[,2]<-1+p_pyrhulla_bed_output[,2]
p_pyrhulla_bed_output
write.table(p_pyrhulla_bed_output, "C:/Users/vik16/OneDrive/Documents/FOGS Summer 2022/Pyrhulla_pyrhulla_PrimerRef_perf.bed",
            row.names=FALSE,col.names=FALSE)
