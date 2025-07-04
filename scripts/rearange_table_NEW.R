# This script is used in TP53 pipeline for table rearangments

library(data.table)
library(outliers)

args = commandArgs(trailingOnly=TRUE)
#setwd('/storage/TP53/data/TP53_20240802-reseq_1/tmp/')
#args = c()

#args[1] = '/storage/TP53/data/TP53_20240802-reseq_1/tmp/'
#args[2] = "TEST"
annot_files = list.files(args[1], pattern = "*normVEP.txt")


print("filter dataset and add information")

# filter dataset and add information
annot_dt_list = lapply(annot_files, function(file){

  # read tabel
  annot_dt = fread(file, header = T)
  annot_dt = cbind(Sample_name = gsub(".norm.merged.annot.normVEP.txt", "", file), annot_dt)


  # variant detected by how many callers?
  annot_dt[set == "vardict-varscan_SNV", VC_DETECTION_COUNT := "2"]
  annot_dt[set == "vardict-varscan_INDEL", VC_DETECTION_COUNT := "2"]
  annot_dt[set == "varscan_SNV", VC_DETECTION_COUNT := "1"]
  annot_dt[set == "varscan_INDEL", VC_DETECTION_COUNT := "1"]
  annot_dt[set == "vardict", VC_DETECTION_COUNT := "1"]

  # filter transcripts
  annot_dt = annot_dt[grep("NM_000546", annot_dt$VEP_Feature),]
  annot_dt$VEP_HGVSp = gsub("%3D", "=", annot_dt$VEP_HGVSp)
  annot_dt = annot_dt[!duplicated(annot_dt$VEP_HGVSc),]

})

    # create big data table from list for all samples
    tmp_dt <- rbindlist(annot_dt_list)
    # remove columns redundant
    rem_col_vector = c("vardict:FREQ","vardict:GQ","vardict:GT",
                      "vardict:PVAL","vardict:RBQ","vardict:RD","vardict:RDF","vardict:RDR",
                      "vardict:SDP","vardict:VD","ID","QUAL","FILTER","AC","ADJAF","ADP","AF","AMPFLAG","AN",
                      "BIAS","CSQ","GDAMP","HET","HICNT","HOM","LSEQ","MQ","MSI","MSILEN","NC","NCAMP",
                      "NM","ODDRATIO","PMEAN","PSTD","QSTD","RSEQ","SAMPLE","SBF","SHIFT3","SN",
                      "TLAMP","TYPE","VARBIAS","VD","VEP_Allele","VEP_BIOTYPE","VEP_DISTANCE",
                      "VEP_Existing_variation","VEP_FLAGS","VEP_Feature_type","VEP_Gene","VEP_HGNC_ID",
                      "VEP_HGVS_OFFSET","VEP_IMPACT","VEP_PolyPhen","VEP_Protein_position","VEP_REFSEQ_MATCH",
                      "VEP_SIFT","VEP_STRAND","VEP_SYMBOL_SOURCE","VEP_cDNA_position","WT","vardict:ABQ","DUPRATE","SPANPAIR","SPLITREAD","SVLEN","SVTYPE")
    tmp_dt[,(rem_col_vector):=NULL]

    # If only vardict call variant, put NA to varscan columns
    for(col_name in grep("varscan", names(tmp_dt),value = T)){
      tmp_dt[set == "vardict",(col_name) := NA]
    }

    # get rid off ":"
    colnames(tmp_dt) = gsub(":", ".", colnames(tmp_dt))
    tmp_dt2 = tmp_dt

    print("join varscan allelic FREQ")
    # join varscan allelic FREQ
    tmp_dt2[,varscan_INDEL.FREQ := gsub("^.$", "F", tmp_dt2$varscan_INDEL.FREQ)]
    tmp_dt2[,varscan.FREQ := tmp_dt2[, .(id=tmp_dt2[["id"]], varscan.FREQ=do.call(paste, c(.SD, sep=""))),
                                 .SDcols= c("varscan_INDEL.FREQ", "varscan_SNV.FREQ")]]
    tmp_dt2[, varscan.FREQ := gsub("%.", "%", tmp_dt2$varscan.FREQ)]
    tmp_dt2[, varscan.FREQ := gsub("F", "", tmp_dt2$varscan.FREQ)]

    # count variant mean across samples
    tmp_dt2[varscan.FREQ != "NANA" & varscan.FREQ != "NoneNone" & varscan.FREQ != ".",final_AF := varscan.FREQ] 
    tmp_dt2[is.na(final_AF), final_AF := as.numeric(as.character(vardict.AF))*100]
    #NoneNone, 0.11%.11% apod
    tmp_dt2[,final_AF := as.numeric(gsub("%.*", "", final_AF))] #tady neco spatne
    tmp_dt2[,variant_mean_perc := round(mean(final_AF), digits = 3),by = VEP_HGVSc]
    tmp_dt2[,variant_median := median(final_AF),by = VEP_HGVSc]
    
    # how often same mutation (AF > 0.05) occur in run
    tmp_dt2[final_AF > 0.05, Occurence_in_samples := .N, by = VEP_HGVSc]

    # mutation description column
#    tmp_dt2[, mutation_description := paste0(gsub("NM_000546.5:","", tmp_dt2$VEP_HGVSc), " ", gsub("NP_000537.3:", "", tmp_dt2$VEP_HGVSp), " ",
#                                                    gsub("\\.",",",as.numeric(as.character(tmp_dt2$vardict.AF)) * 100), "%")]
#    tmp_dt2[varscan.FREQ != "NANA", mutation_description := paste0(gsub("NM_000546.5:","", tmp_dt2$VEP_HGVSc), " ", gsub("NP_000537.3:", "", tmp_dt2$VEP_HGVSp), " ",
#                                                    gsub("\\.", ",", tmp_dt2$varscan.FREQ))]
    for(p in 1:length(tmp_dt2$Sample_name)){

      if(tmp_dt2$varscan.FREQ[p] == "NANA"){
        tmp_dt2$mutation_description[p] = paste0(gsub("NM_000546.5:","", tmp_dt2$VEP_HGVSc[p]), " ", gsub("NP_000537.3:", "", tmp_dt2$VEP_HGVSp[p]), " ",
                                                    gsub("\\.",",",as.numeric(as.character(tmp_dt2$vardict.AF[p])) * 100), "%")
      }else{
        tmp_dt2$mutation_description[p] = paste0(gsub("NM_000546.5:","", tmp_dt2$VEP_HGVSc[p]), " ", gsub("NP_000537.3:", "", tmp_dt2$VEP_HGVSp[p]), " ",
                                                    gsub("\\.", ",", tmp_dt2$varscan.FREQ[p]))
      }
    }

    print("count strand bias")
    # count strand bias
    for(r in 1:length(tmp_dt2$Sample_name)){
      # all Varscan variants
      if(tmp_dt2$set[r] == paste0("varscan_SNV") || tmp_dt2$set[r] == paste0("varscan_INDEL") || tmp_dt2$set[r] == paste0("vardict-varscan_INDEL") ||
         tmp_dt2$set[r] == paste0("vardict-varscan_SNV") || tmp_dt2$set[r] == paste0("filterInvardict-varscan_SNV")){
        
        # all Varscan SNV
        if(tmp_dt2$set[r] == paste0("varscan_SNV") || tmp_dt2$set[r] == paste0("vardict-varscan_SNV") || tmp_dt2$set[r] == paste0("filterInvardict-varscan_SNV")){
          #print("SNV")
          citatel <- as.numeric(as.character(tmp_dt2$varscan_SNV.ADF[r]))/(as.numeric(as.character(tmp_dt2$varscan_SNV.ADF[r]))+as.numeric(as.character(tmp_dt2$varscan_SNV.RDF[r])))
          jmenovatel <- as.numeric(as.character(tmp_dt2$varscan_SNV.ADR[r]))/(as.numeric(as.character(tmp_dt2$varscan_SNV.ADR[r]))+as.numeric(as.character(tmp_dt2$varscan_SNV.RDR[r])))
          tmp_dt2$varscan.SBIAS[r] <- round(citatel/jmenovatel, digits = 2)
        }else{
          # all Varscan INDEL
          #print("INDEL")
          citatel <- as.numeric(as.character(tmp_dt2$varscan_INDEL.ADF[r]))/(as.numeric(as.character(tmp_dt2$varscan_INDEL.ADF[r]))+as.numeric(as.character(tmp_dt2$varscan_INDEL.RDF[r])))
          jmenovatel <- as.numeric(as.character(tmp_dt2$varscan_INDEL.ADR[r]))/(as.numeric(as.character(tmp_dt2$varscan_INDEL.ADR[r]))+as.numeric(as.character(tmp_dt2$varscan_INDEL.RDR[r])))
          tmp_dt2$varscan.SBIAS[r] <- round(citatel/jmenovatel, digits = 2)
        }
      }else{
        tmp_dt2$varscan.SBIAS[r] <- "NA"
      }
    }
    tmp_dt2[,varscan.SBIAS := gsub("Inf", "10000",tmp_dt2$varscan.SBIAS)]
    tmp_dt2[,varscan.SBIAS := gsub("NA", "4.44", tmp_dt2$varscan.SBIAS)]
    
    tmp_dt2$varscan.SBIAS =  as.numeric(tmp_dt2$varscan.SBIAS)
    

    ######################## evaluate background variants ######################################################
    print("evaluate background variants")
    tmp_dt2[is.na(tmp_dt2$final_AF)]
    #nektere final_AF jsou NULL
    tmp_dt2[,outlier_p_value := .(grubbs.test(final_AF)$p.value), by = VEP_HGVSc]
    tmp_dt2[outlier_p_value < 0.05, max_value := .(max(final_AF)), by = VEP_HGVSc]
    
    
    tmp_dt2[,vyskyt := .(round(Occurence_in_samples/length(unique(Sample_name)), digits = 3))]
    
    tmp_dt2[is.na(max_value), sd275 := round(2.75*sd(final_AF), digits = 3), by = VEP_HGVSc]
    tmp_dt2[!is.na(max_value), sd275 := round(2.75*sd(final_AF[!final_AF %in% max_value[1]]), digits = 3), by = VEP_HGVSc]
    
    tmp_dt2[is.na(max_value), vrmean := round(mean(final_AF), digits = 3), by = VEP_HGVSc]
    tmp_dt2[!is.na(max_value), vrmean := round(mean(final_AF[!final_AF %in% max_value[1]]), digits = 3), by = VEP_HGVSc]
    
    tmp_dt2[vyskyt > 0.75 & final_AF < vrmean + sd275, background := "YES_sarka_rule"]

    # Extract rows with more than 1 ALT
    matching_rows <- grepl(",", tmp_dt2$ALT)
    tmp_dt2$background[matching_rows] <- "YES_sarka_rule"
    
    tmp_dt2[varscan.SBIAS < 0.2 | varscan.SBIAS > 5, background := "YES_strandbias"]
    
    
    # Set variants detected only by vardict to NA (previously falsely marked as YES_strandbias)
    tmp_dt2[set == "vardict", background := NA]
    
    # Step 1: Identify rows where VEP_HGVSc contains "c.215C>G"
    matching_rows <- grepl("c.215C>G", tmp_dt2$VEP_HGVSc)
    
    # Step 2: Add rows where VEP_HGVSc contains "c.108G>A"
    matching_rows <- matching_rows | grepl("c.108G>A", tmp_dt2$VEP_HGVSc)
    
    # Step 3: Add rows where VEP_HGVSc contains "c.639A>G"
    matching_rows <- matching_rows | grepl("c.639A>G", tmp_dt2$VEP_HGVSc)
    
    # Rename the "background" column values to "polymorfismus" for these rows
    tmp_dt2$background[matching_rows] <- "polymorfismus"
  
    
    # Fix missing QUALS for non-Vardict variants
    # Step 1: Extract the 'VarScanType' directly in the data.table
    tmp_dt2[, VarScanType := tstrsplit(set, "_")[[2]]]
    
    # Step 2: Use a vectorized approach to update 'QUAL' only if a condition is met
    tmp_dt2[, QUAL := fifelse(VarScanType == "INDEL", varscan_INDEL.ABQ, 
                              fifelse(VarScanType == "SNV", varscan_SNV.ABQ, QUAL))]
    
    
    ##########################################################################################################
    
    tmp_dt2[,varscan.SBIAS := gsub("4.44", "NA", tmp_dt2$varscan.SBIAS)]
    print("reformat table")

    # reformat table
    tmp_dt2[,VEP_INTRON := gsub("/10", "", tmp_dt2$VEP_INTRON)]
    tmp_dt2[,VEP_EXON  := gsub("/11", "", tmp_dt2$VEP_EXON)]
    tmp_dt2[,vardict.AD := gsub(",",";", tmp_dt2$vardict.AD)]
    tmp_dt2[,vardict.AF := paste0(gsub("\\.",",",as.numeric(as.character(tmp_dt2$vardict.AF)) * 100))]
    colnames(tmp_dt2)[colnames(tmp_dt2) == 'vardict.AF'] <- 'vardict.AF_percentage'
    tmp_dt2[,vardict.ALD := gsub(",",";",tmp_dt2$vardict.ALD)]
    # fix numeric vals
    tmp_dt2$HIAF = as.numeric(tmp_dt2$HIAF)
    tmp_dt2$varscan.FREQ = as.numeric(sub("%", "", tmp_dt2$varscan.FREQ))
    tmp_dt2$varscan_INDEL.PVAL = as.numeric(tmp_dt2$varscan_INDEL.PVAL)
    tmp_dt2$varscan_SNV.PVAL = as.numeric(tmp_dt2$varscan_SNV.PVAL)
    tmp_dt2$variant_mean_perc = as.numeric(tmp_dt2$variant_mean_perc)
    tmp_dt2$variant_median = as.numeric(tmp_dt2$variant_median)
    tmp_dt2$varscan.SBIAS = as.numeric(tmp_dt2$varscan.SBIAS)
    
    # tmp_dt2[,HIAF := gsub("\\.",",", tmp_dt2$HIAF)]
    # tmp_dt2[,varscan.FREQ := gsub("\\.", ",", tmp_dt2$varscan.FREQ)]
    # tmp_dt2[,varscan_INDEL.PVAL := gsub("\\.", ",", tmp_dt2$varscan_INDEL.PVAL)]
    # tmp_dt2[,varscan_SNV.PVAL := gsub("\\.", ",", tmp_dt2$varscan_SNV.PVAL)]
    # tmp_dt2[,variant_mean_perc := gsub("\\.", ",", tmp_dt2$variant_mean_perc)]
    # tmp_dt2[,variant_median := gsub("\\.", ",", tmp_dt2$variant_median)]
    # tmp_dt2[,varscan.SBIAS := gsub("\\.", ",", tmp_dt2$varscan.SBIAS)]

# create and save final tables
DIR_file_name <- paste0(args[2],".allsamples.merged.anot.txt")
write.table(tmp_dt2,file=DIR_file_name, sep="\t", quote = F, row.names = F)

final_dt = tmp_dt2
for(sample in unique(final_dt$Sample_name)){
  # subset samples
  # print(final_dt$Sample_name)
  final_dt2 = final_dt[grep(sample, final_dt$Sample_name),]
  # remove from each file
  remove_columns = c("DP", "HIAF", "HICOV", "QUAL.1", "REFBIAS", "vardict.DP", "vardict.ADF", "vardict.ADR", "varscan_INDEL.ABQ", "varscan_INDEL.AD",
                     "varscan_INDEL.AF", "varscan_INDEL.ALD",
                     "varscan_INDEL.GQ", "varscan_INDEL.GT", "varscan_INDEL.PVAL", "varscan_INDEL.RBQ",
                     "varscan_INDEL.RD", "varscan_INDEL.SDP",
                     "varscan_INDEL.VD", "varscan_INDEL.FREQ", "varscan_SNV.FREQ", "varscan_SNV.ABQ", "varscan_SNV.AD",
                     "varscan_SNV.AF", "varscan_SNV.ALD", "varscan_SNV.GQ", "varscan_SNV.GT", "varscan_SNV.PVAL",
                     "varscan_SNV.RBQ", "varscan_SNV.RD", "varscan_SNV.SDP",
                     "varscan_SNV.VD", "VC_DETECTION_COUNT", "final_AF")
  final_dt2 = final_dt2[,(remove_columns):=NULL]
  # save file
  print(paste0(sample, ".sample.merged.anot.txt"))
  #print(head(final_dt2))
  write.table(final_dt2, file=paste0(sample, ".sample.merged.anot.txt"), sep="\t", quote = F, row.names = F)
}

