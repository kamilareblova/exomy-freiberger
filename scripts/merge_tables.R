library(data.table)
#library(tidyr)
#library(tictoc)

# !!! uncomment
#tic("count strand bias")

args = commandArgs(trailingOnly=TRUE)

# arguments temporary

# Path to inputs
print(args[1])
# RUN ID
print(args[2])

merge_l = list()
for(i in list.files(args[1], pattern = ".sample.merged.anot.txt")){

	id = gsub(".sample.merged.anot.txt", "", i)
	dt = fread(paste0(args[1], i))	
	dt[,sampleID := id]
	merge_l[[i]] = dt

}

final_dt = rbindlist(merge_l)
write.table(final_dt, file = paste0(args[1], "/", args[2], ".allsamples.merged.anot.txt"), sep = "\t", quote = F, row.names = F)
