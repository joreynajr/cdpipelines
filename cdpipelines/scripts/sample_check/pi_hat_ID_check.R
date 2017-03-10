#!/usr/bin/env /software/R-3.2.2-cardips/bin/Rscript

library(ggplot2)
library(stringr)
suppressPackageStartupMessages(library(RMySQL))
source("/frazer01/home/matteo/my_software/cardips_functions.R")
source("/frazer01/home/matteo/my_software/ase_pipeline/ase_pipeline_functions.R") 
source("/nas/Personal folders/margaret/interintra_twins/scripts/twin_functions.R")
source("/nas/Personal folders/margaret/interintra_twins/scripts/ase_pipeline_functions_md.R")

# functions #
checkTissue = function(col1, col2){
    
    if (col1 == 'iPSC' | col1 == "CM"){
        
        return(col1)
    }
    else if (col2 == 'iPSC' | col2 == "CM"){
        
        return(col2)
    }
    
    else{
        return("not_iPSC_or_CM_in_db")
    }
}

checkMatch = function(one, two){
   
    x = as.integer(substr(one, 2, nchar(one)))
    y = as.integer(substr(two, 2, nchar(two)))
    
    if (x == y){
        
        return(1)
    }
    
    else if (x == (y + 1)){
        
        return(1)
    }
    
    else if (x == (y - 1)){
        
        return(1)
    }
    
    else {
        
        return(0)
    }
    
    
}

## 

query_rna = "SELECT * FROM data_rnas;"
db_table_rna_data = getTableFromDbMD(query_rna)
db_table_rna_data$id = unlist(lapply(db_table_rna_data$id, addDashUUIDMD))
filt_all_rna <- data.frame(db_table_rna_data$id, str_split_fixed(db_table_rna_data$name, "_", 6)[,1])
colnames(filt_all_rna) <- c('uuid', 'indv')

query_wgs = "SELECT * FROM data_wgs;"
db_table_wgs = getTableFromDbMD(query_wgs)
db_table_wgs$id = unlist(lapply(db_table_wgs$id, addDashUUIDMD))
filt_db_table_wgs <- data.frame(db_table_wgs$id, str_split_fixed(db_table_wgs$name, "_", 6)[,1])
colnames(filt_db_table_wgs) <- c('uuid', 'indv')

tissue_df <- data.frame(db_table_rna_data$id, str_split_fixed(db_table_rna_data$name, "_", 6)[,2], str_split_fixed(db_table_rna_data$name, "_", 6)[,3])

colnames(tissue_df) <- c("uuid", "one", "two")
tissue_df$tissue <- tissue_df$two
tissue_df$one <- sub("^$", "N", tissue_df$one)
tissue_df$two <- sub("^$", "N", tissue_df$two)
tissue_df$tissue <- sub("^$", "N", tissue_df$tissue)


for(row in seq(1, nrow(tissue_df))){
    tissue_df[row, 4] = checkTissue(as.character(tissue_df[row, 2]), as.character(tissue_df[row, 3]))
}

tissue_df$one <- NULL
tissue_df$two <- NULL


all_plink_files <- Sys.glob("/projects/CARDIPS/pipeline/RNAseq/sample/*/qc/*_sample_swap/plink/merged.genome.sample")
l <- length(all_plink_files)
master_match <- data.frame(seq(1, l), seq(1, l), seq(1, l), seq(1, l), seq(1, l), seq(1, l), row.names = seq(1, l))
colnames(master_match) <- c("sample_tissue", "sample_dataid", "sample_subject", "WGSmatch_dataid", "WGSmatch_subject", "match_status")

counter = 0

for (plink in all_plink_files){

    uuid = str_split_fixed(plink, "/", 11)[,7]
    path_to_out = paste(str_split_fixed(plink, "/", 11)[,1:9], collapse = "/")
    counter = counter + 1

    data <- read.table(plink, header = T)
    filt_data <- data[data$FID1 == uuid | data$FID2 == uuid, ]
    maxRow <- (filt_data[filt_data$PI_HAT == max(filt_data$PI_HAT),])

    sample <- (as.character(maxRow$IID1)[1]) #max might match more than one sample, because MO twins and clones are included in analysis. Take first one.
    ref <- (as.character(maxRow$IID2))[1] #max might match more than one sample, because MO twins and clones are included in analysis. Take first one.
    
    # hist(filt_data$PI_HAT, xlab = "PI_HAT", main = uuid)

    # ## filt_data uuid might occur in IID1 or IID2 columns, so I need to check both columns.
    s <- (filt_all_rna[filt_all_rna$uuid == sample,])
    r <- (filt_db_table_wgs[filt_db_table_wgs$uuid == ref,])

    s_swap <- (filt_all_rna[filt_all_rna$uuid == ref,])
    r_swap <- (filt_db_table_wgs[filt_db_table_wgs$uuid == sample,])

    tissue = "iPSC"
    ## If it did not occur in in IID1, it must have occured in IID2, so IID1 is the WGS reference data uuid
    if (length(as.character(s$indv)) == 0){

        match = checkMatch(as.character(s_swap$indv), as.character(r_swap$indv))
        master_match[counter, 1] = tissue_df[tissue_df$uuid == uuid, 2]
        master_match[counter, 2] = ref
        master_match[counter, 3] = as.character(s_swap$indv)
        master_match[counter, 4] = sample
        master_match[counter, 5] = as.character(r_swap$indv)
        master_match[counter, 6] = as.character(match)
        

    }

    ## It occured in IID1, so IID2 is the WGS reference data uuid
    if (length(as.character(s$indv)) == 1){

        match = checkMatch(as.character(s$indv), as.character(r$indv))
        
        master_match[counter, 1] = tissue_df[tissue_df$uuid == uuid, 2]
        master_match[counter, 2] = sample
        master_match[counter, 3] = as.character(s$indv)
        master_match[counter, 4] = ref
        master_match[counter, 5] = as.character(r$indv)
        master_match[counter, 6] = as.character(match)
        
    }
}

write.csv(master_match, file="/frazer01/home/joreyna/repos/cdpipelines/cdpipelines/scripts/sample_check/master_rna_sample_identity.csv", quote = F, row.names = F)
