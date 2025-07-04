library(data.table)
args = commandArgs(trailingOnly=TRUE)
path<-args[1]
#path<-"/mnt/nfs/home/357156/999993-Bioda/projects/TP53/20180625AndreaMareckova/coverage"

final_list<-list()
cov_list_tables = list()

for(i in list.files(path, pattern = "*PBcoverage.txt")){
  
  # table of per base coverage for each sample (to paste together to one table)
  file <- read.table(paste0(path,"/",i), sep="\t")
  colnames(file) <- c("Chrom", "Exon", "K_NICEMU", "Start", "End", "K_NICEMU2", "STRAND", "K_NICEMU3", "Base_number", "Coverage")
  file <- file[,-c(3,6,8)]
  file$sample <- gsub(".PBcoverage.txt", "", i)
  final_list[[i]] <- file
  
  # table of stats per each sample
  file_name_stat <- gsub(".PBcoverage.txt", "", i)
  
  exon_stat <- list()
  for(j in levels(file$Exon)){
    # only particular exon in particular sample
    tmp<-file[grep(j,file$Exon),]
    
    tmp_final <- data.frame(Chrom = tmp$Chrom[1], Exon = tmp$Exon[1], Start = tmp$Start[1], End = tmp$End[1],
                            Coverage_max = max(tmp$Coverage), Coverage_min = min(tmp$Coverage),
                            Coverage_mean = round(mean(tmp$Coverage), digit = 2), Coverage_median = median(tmp$Coverage))
    exon_stat[[j]]<-tmp_final
  }
  
  exon_sample_stat <- rbindlist(exon_stat)
  write.table(exon_sample_stat, file = paste0(file_name_stat,".perexon_stat.txt"), quote = F, row.names = F, sep = "\t")

	exon_sample_stat[,sampleID := file_name_stat]
	cov_list_tables[[i]] = exon_sample_stat}

final_df<-rbindlist(final_list)
write.table(final_df, file=paste0(args[2], "_tp53_gene_coverage.txt"),
            sep="\t", quote = F, row.names = F)

final_dt_stats = rbindlist(cov_list_tables)
write.table(final_dt_stats, file = paste0(args[2], ".allsamples.merged.coverage.txt"), sep = "\t", quote = F, row.names = F)

