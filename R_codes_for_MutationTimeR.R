### R codes for running MutationTimeR using 
rm(list = ls())

# Install dependency packages and load in data
library(tidyverse)
library(VariantAnnotation)
library(MutationTimeR)

load("ITH2_data.Rdata")
ITH2.table <- read.csv("Table.S4.csv")
ITH2.hatchet <- read.csv("data-SCNA-segment-HATCHet.csv") %>% 
  subset(ID %in% ITH2.clinical$Patient_ID)

ITH2.mut.cat$Driver.gene <- ifelse(ITH2.mut.cat$Gene %in% ext.data.driver.list, 
                                   "Driver", "Non-driver")


# Main function of mutationtimer, estimating clonality of a somatic mutation based on 
# the VAF, tumor purity, and A/B copy number of that allele

mut.timing <- function(ID){
  # prepare input data needed for running MutationTimeR
  mut.cat <- subset(ITH2.mut.cat, Sector_ID == ID & Chromosome %in% 1:22) %>% 
    mutate(Alt_c = round(Depth_tumor_sample*Variant_allele_fraction, 0)) %>% 
    mutate(Ref_c = Depth_tumor_sample - Alt_c) %>% 
    mutate(mutation = str_glue("{Chromosome}:{Start_position}_{Ref}/{Alt}"))
  
  PT <- unique(mut.cat$Patient_ID)
  purity <- ITH2.table$Tumor_purity[ITH2.table$Sector_WES_ID == ID] %>% as.numeric()
  ploidy <- ITH2.table$Tumor_ploidy[ITH2.table$Sector_WES_ID == ID] %>% as.numeric()
  seg <- subset(ITH2.hatchet, ID == PT) %>% 
    mutate(CHR = sub("chr","", CHR))
  WGD.status <- ifelse(unique(seg$BG.ploidy)==4, TRUE, FALSE)
  
  cat(paste("Timing gene mutation for", PT, ID, "\n"))
  cat(paste("Tumor purity:", purity, "\n"))
  cat(paste("Genome doubling status:", WGD.status, "\n"))
  
  # create a CN segment dataframe at the arm level as an input for MutationTimeR
  tmp.path <- paste0("./",ID,".vcf")
  CNV.region <- read.delim("data-GRCh38p7-region.tsv")
  CNV.region$chr <- sub("chr","",x = CNV.region$chr)
  seg.arm <- seg
  seg.arm$CHR <- sapply(1:nrow(seg.arm), function(x){
    arm <- ifelse(seg.arm$START[x] < CNV.region$centromere_start[CNV.region$chr == seg.arm$CHR[x]],
                  "", ".5")
    return(paste0(seg.arm$CHR[x],arm))
  })
  
  tmp <- lapply(unique(seg.arm$CHR), function(x){
    df <- subset(seg.arm, CHR == x) 
    seg.chr <- data.frame(CHR = x, start.pos = min(df$START), end.pos = max(df$END),
                          major_CN = round(sum(df$CNa*(df$END-df$START))/sum(df$END-df$START, digits = 0)),
                          minor_CN = round(sum(df$CNb*(df$END-df$START))/sum(df$END-df$START, digits = 0)))
    return(seg.chr)
  })
  tmp <- do.call("rbind", tmp)
  tmp$CHR <- sub("\\.5","",tmp$CHR) %>% sub(pattern = "chr", replacement = "")
  
  gr.arm <- GRanges(Rle(c(tmp$CHR)), 
                    IRanges(start = tmp$start.pos, end = tmp$end.pos),
                    major_cn=as.integer(tmp$major_CN), minor_cn=as.integer(tmp$minor_CN), 
                    clonal_frequency=purity)
  
  # create a vcf file and load in as input for MutationTimeR
  vcf.tmp <- data.frame(Chr=mut.cat$Chromosome, Start=mut.cat$Start_position, ID = ".", 
                        Ref=mut.cat$Ref, Alt=mut.cat$Alt, Qual = ".", FILTER = "PASSED",
                        info = paste0("t_alt_count=", mut.cat$Alt_c,";t_ref_count=", mut.cat$Ref_c),
                        Format = "AD:DP", 
                        Sample = paste0(".,",mut.cat$Alt_c,":",mut.cat$Depth_tumor_sample))
  
  header <- readLines("./MutationTimeR.vcf") %>% grep(pattern = "#", value = T)
  cat(header, file = tmp.path, sep = "\n")
  write.table(vcf.tmp[,1:10], file = tmp.path, append = T, quote = F, sep = "\t", row.names = F, col.names = F)
  
  vcf <- VariantAnnotation::readVcf(tmp.path) 
  unlink(tmp.path)
  
  # Run MutationTimeR with n.boot = 1000
  clusters <- data.frame(cluster = 1:2, proportion = purity*c(1, 0.5), n_ssms=c(90,10))
  
  mt.arm <- mutationTime(vcf = vcf, cn = gr.arm, isWgd = WGD.status, n.boot=1000, 
                         clusters = clusters, purity = purity)
  
  vcf <- addMutTime(vcf, mt.arm$V)
  vcf.timing <- info(vcf) %>% data.frame()
  vcf.timing <- mutate(vcf.timing, mutation = rownames(vcf.timing), ID = ID, 
                       WGD = WGD.status, purity = purity, ploidy = ploidy) %>% 
    left_join(y = dplyr::select(mut.cat, mutation, Gene, Amino_acid_change, Exonic_function,
                                Variant_allele_fraction, Depth_tumor_sample),
              by = "mutation") %>% 
    arrange(CLS, desc(Variant_allele_fraction))
  
  return(vcf.timing)
  
}
tmp <- mut.timing(ID = "CHL664")

Mut.timing <- lapply(unique(ITH2.mut.cat$Sector_ID), mut.timing)
Mut.timing <- do.call("rbind", Mut.timing)
Mut.timing <- left_join(Mut.timing, dplyr::select(ITH2.table, Patient_ID, Sector_WES_ID), 
                        by = c("ID"="Sector_WES_ID"))

