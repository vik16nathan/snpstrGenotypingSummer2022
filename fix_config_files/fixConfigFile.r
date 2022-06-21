library(readr)

#Take the primer as an input
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=1) {
  stop("Must supply primer only!", call.=FALSE)
} 

primer <- args[1]

#Load in the .config file
old_config_file <- as.data.frame(read_tsv(paste(
    "./config/Pyrhulla_pyrhulla_Primer",primer,
    "_multiple_flanks.config",sep="")))

fixed_flanks_file <- as.data.frame(read_tsv(paste(
    "Primer",primer,"_corrected_flanks.tsv",sep="")))


for(row in c(1:nrow(old_config_file))) {

    if(old_config_file[row,"5'Flank"] %in% fixed_flanks_file[,"Incorrect Flank"]) {
        old_config_file[row,"5'Flank"] <- fixed_flanks_file[
            which(fixed_flanks_file[,"Incorrect Flank"]==old_config_file[row,"5'Flank"]),"Corrected Flank"]
    }

    if(old_config_file[row,"3'Flank"] %in% fixed_flanks_file[,"Incorrect Flank"]) {
        old_config_file[row,"3'Flank"] <- fixed_flanks_file[
            which(fixed_flanks_file[,"Incorrect Flank"]==old_config_file[row,"3'Flank"]),"Corrected Flank"]
    }
    
}

old_config_file <- unique(old_config_file)
write.table(old_config_file,paste(
    "./config/Pyrhulla_pyrhulla_Primer",primer,
    "_multiple_flanks.config",sep=""),sep="\t",quote=FALSE,row.names=FALSE)