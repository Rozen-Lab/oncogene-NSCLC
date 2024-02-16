### R codes for generating input files for pyclone-iv
rm(list = ls())
library(tidyverse)

load("ITH2_data.Rdata")
ITH2.table <- read.csv("Table.S4.csv") %>% 
  subset(Tissue_type == "Lung tumor" & Tumor_purity != "< 0.1")
mut.timing <- read.csv("data-mutation-timing-clonality.csv")

lapply(unique(ITH2.table$Patient_ID), function(x){
  output.path <- paste0("./pyclone/input/", x, ".tsv")
  cat(paste0("writing table for ", x, "\n"))

  ID.df <- subset(ITH2.table, Patient_ID == x)$Sector_WES_ID
  mut.cat <- subset(mut.timing, Sector.ID %in% ID.df) 
  df <- dplyr::select(mut.cat, mutation_id = Gene.mut, sample_id = Sector.ID, 
                      ref_counts = Tumor_ref_count, alt_counts = Tumor_alt_count, 
                      normal_cn = 2, major_cn = MajCN, minor_cn = MinCN, 
                      tumour_content = Tumor.purity)
  df <- df[!is.na(df$major_cn) & !is.na(df$minor_cn), ]
  
  write.table(df, file = output.path, sep = "\t", col.names = T, row.names = F, quote = F)
})

