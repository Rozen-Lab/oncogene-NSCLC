### R codes for oncoplot & stats testing for clonality (Fig.1A)
# Install dependency packages
rm(list = ls())
library(dplyr)
library(tidyr)
library(magrittr)
library(stringr)
load("ITH2_data.Rdata")

# calculate tumor mutational burden per sector (defined by the number of 
# unique mutations) and intra-tumor heterogeneity per sector (defined by the 
# proportion of branch to total mutations)
ITH2.mut.cat <- mutate(ITH2.mut.cat, Gene.mut = str_glue("{Gene}_{Amino_acid_change}"))

tmp <- count(ITH2.mut.cat, Sector_ID, Clonality) %>% 
  spread(Clonality,n,fill=0) %>% 
  mutate(TMB.sector = Branch + Truncal, ITH.sector = Branch / (Branch + Truncal)) %>% 
  mutate(Branch.prop = Branch / TMB.sector, Truncal.prop = Truncal / TMB.sector)

ITH2.table <- subset(ITH2.table, Tissue_type == "Tumor") 

ITH2.table <- left_join(x = dplyr::select(ITH2.table, Group, ends_with("ID")), y = tmp, 
                        by = c("Sector_WES_ID" = "Sector_ID")) %>% 
  left_join(y = dplyr::select(ITH2.clinical, Patient_ID, Gender, Tumor_histology, 
                              Smoking_pk_yr, Oncogene_mut), by = "Patient_ID")

tmp2 <- group_by(ITH2.table, Patient_ID) %>% 
  summarize(TMB.max = max(TMB.sector), TMB.tumor = mean(TMB.sector))

ITH2.table <- left_join(x = ITH2.table, y = tmp2, by = "Patient_ID")

# calculate the prevalence of selected driver mutation in the cohort
tmp <- subset(ITH2.mut.cat, Driver_gene == "Driver") %>% 
  mutate(CHECK = str_glue("{Patient_ID}-{Gene}")) %>% 
  subset(!duplicated(CHECK)) %>% 
  count(Patient_ID, Gene) %>% 
  spread(Gene,n,fill=0)

gene.freq <- data.frame(gene = colnames(tmp)[-1], mutation = colSums(tmp[,-1]), N = nrow(tmp))
gene.freq$rate <- round(gene.freq$mutation*100/gene.freq$N, digits = 1) %>% format(nsmall = 1)
gene.freq <- arrange(gene.freq, desc(mutation))

# get the status of driver mutations
gene.highlight <- c(
  "EGFR","MET","ERBB2","KRAS","TP53","CSMD3","MUC16","RB1","PTPRD","EPHA3","CNBD1","CTNND2",
  "RBM10","BIRC6","CTNNB1","FKBP9","MTOR","KEAP1","CASP9","APC","STK11")

ITH2.table$TP53.status <- 
  ifelse(ITH2.table$Patient_ID %in% subset(ITH2.mut.cat, Gene == "TP53")$Patient_ID, 
         "Mutant", "Wild type")

ITH2.table$Group <- factor(
  ITH2.table$Group, 
  levels = c("Oncogene-driven non-smoking","Oncogene-driven smoking","Typical smoking"))

ITH2.table$Oncogene_mut <- factor(
  ITH2.table$Oncogene_mut, 
  levels = c("EGFR","MET","ERBB2","ALK","KRAS","Wild type"))

ITH2.table$RNA.subtype <- ifelse(is.na(ITH2.table$Sector_RNA_ID), "No RNA sample",
                                 ifelse(ITH2.table$Sector_RNA_ID %in% cluster.TRU, 
                                        "TRU", "Non-TRU"))

ITH2.table <- arrange(ITH2.table, Group, Oncogene_mut, TP53.status, 
                      desc(TMB.tumor), desc(TMB.max), desc(TMB.sector))

### add gene fusion from STAR-Fusion output to the driver gene list
fusion.highlight <- c("EML4--ALK", "KLC1--ALK", "PARG--BMS1")

fusion.label <- sapply(ITH2.table$Sector_WES_ID, function(x){
  if(is.na(ITH2.table$Sector_RNA_ID[ITH2.table$Sector_WES_ID==x])==TRUE){
    Status <- rep("No RNA sample", length(fusion.highlight))
  }
  else if(is.na(ITH2.table$Sector_RNA_ID[ITH2.table$Sector_WES_ID==x])==FALSE){
    tmp <- subset(ITH2.gene.fusion, Sector_ID == x & Fusion_gene %in% fusion.highlight)
    Status <- sapply(fusion.highlight, function(x){
      ifelse(x %in% tmp$Fusion_gene, subset(tmp, Fusion_gene == x)$Clonality[1], "Wild type")
    })
  } 
  return(Status)
})
colnames(fusion.label) <- ITH2.table$Sector_WES_ID
rownames(fusion.label) <- c("EML4-ALK", "KLC1-ALK", "PARG-BMS1")

# prepare space value of the oncoplot
space.brp <- sapply(unique(ITH2.table$Group), function(x){
  tmp <- subset(ITH2.table, Group == x)$Patient_ID %>% table()
  tmp <- tmp[unique(subset(ITH2.table, Group == x)$Patient_ID)]
  brp <- sapply(tmp, function(x){c(0.5,rep(0,(x-1)))}, simplify = T) %>% unlist()
  brp[1] <- 4
  return(brp)
}) %>% unlist()

# prepare the color key of the oncoplot
color.fusion <- c("No RNA sample" = "white", "Wild type" = "grey95",
                  "Branch" = "lightskyblue", "Truncal" = "blue4")

barplot.mut <- function(gene, plot.mut.type = FALSE){
  mut.n <- subset(ITH2.mut.cat, Gene == gene) %>% count(Sector_ID, Gene)
  clone.n <- max(mut.n$n)
  
  color.mut.v1 <- c("Truncal"="blue4", "Branch"="lightskyblue", 
                    "Mixed"="lightskyblue" ,"Wild type"="grey95")
  color.mut.v2 <- c("Truncal"="blue4", "Branch"="lightskyblue", 
                    "Mixed"="blue4" , "Wild type"="grey95")

  gene.mut <- sapply(ITH2.table$Sector_WES_ID, function(x){
    tmp <- subset(ITH2.mut.cat, Sector_ID == x & Gene == gene)
    if(nrow(tmp)==0){
      return("Wild type")
    } else if(nrow(tmp)>0 & n_distinct(tmp$Clonality)==1){
      return(unique(tmp$Clonality))
    } else if(nrow(tmp)>0 & n_distinct(tmp$Clonality)>1){
      return("Mixed")
    }
  })
  barplot(rep(1,nrow(ITH2.table)), border = NA, axes = FALSE, space = space.brp, 
          col = color.mut.v1[gene.mut])
  barplot(rep(0.5,nrow(ITH2.table)), border = NA, axes = FALSE, space = space.brp, 
          col = color.mut.v2[gene.mut], add = T)
  mtext(gene, side=2, line = 3, cex = 1, las = 2, outer = F, col = 'black')
  mtext(paste0(gene.freq[gene,]$rate, " %"), 
        side=2, line = -2, cex = 1, las = 2, outer = FALSE, col = 'black')
}

## oncoplot 
layout(matrix(c(rep(1,10), rep(2,4),3:30,rep(31,3)), ncol = 1))
par(mar = c(0.2, 0.2, 0.1, 0.2), oma = c(1, 10, 2, 0.5))

#1 plot TMB bar
barplot(as.matrix(t(ITH2.table[,c("Truncal","Branch")])), axes = FALSE, 
        col = c("blue4","steelblue1"), space = space.brp, border = NA, ylim = c(0,800))
axis(side = 2, at = seq(0,800,200), labels = seq(0,800,200), line = -1, las = 1)
mtext("Tumor\nmutation\nburden", side=2, line = 4, cex = 1, las = 2, outer = F, col = 'black')

#2 plot ITH bar
barplot(as.matrix(t(ITH2.table[,c("Branch.prop","Truncal.prop")])), border = NA, axes = FALSE, 
        col = c("gray10","gray80"), space = space.brp, ylim = c(0,1.1))
axis(side = 2, at = c(0,0.5,1), labels = c(0,0.5,1), line = -1, las = 1)
mtext("Intratumoral\nheterogeneity\n(ITH)", side=2, line = 4, cex = 1, las = 2, outer = F, col = 'black')

#3 plot Gender bar
barplot(rep(1,nrow(ITH2.table)), border = NA, axes = FALSE, space = space.brp,
        col = ifelse(ITH2.table$Gender == "Male", "turquoise4", "orchid1"))
mtext("Gender", side=2, line = -2, cex = 1, las = 2, outer = FALSE, col = 'black')

#4 plot smoking bar
barplot(rep(1,nrow(ITH2.table)), border = NA, axes = FALSE, space = space.brp, 
        col = case_when(ITH2.table$Smoking_pk_yr == 0 ~ "grey95",
                        ITH2.table$Smoking_pk_yr > 0   & ITH2.table$Smoking_pk_yr < 20 ~ "grey65",
                        ITH2.table$Smoking_pk_yr >= 20 & ITH2.table$Smoking_pk_yr < 40 ~ "grey35",
                        ITH2.table$Smoking_pk_yr >= 40 ~ "grey5"))
mtext("Smoking", side=2, line = -2, cex = 1, las = 2, outer = FALSE, col = 'black')

#5 plot histology bar
barplot(rep(1,nrow(ITH2.table)), border = NA, axes = FALSE, space = space.brp, 
        col = ifelse(ITH2.table$Tumor_histology == "Adenocarcinoma", "lightblue2", "seagreen4"))
mtext("Histology", side=2, line = -2, cex = 1, las = 2, outer = FALSE, col = 'black')

#6 plot RNA expression subtype
barplot(rep(1,nrow(ITH2.table)), axes = FALSE, space = space.brp, 
        border = ifelse(ITH2.table$RNA.subtype == "No RNA sample", "gray", NA), 
        col = ifelse(ITH2.table$RNA.subtype == "Non-TRU", "steelblue1",
                     ifelse(ITH2.table$RNA.subtype == "TRU", "violetred1", "white")))
mtext("RNA expression subtype", side=2, line = -2, cex = 1, las = 2, outer = FALSE, col = 'black')

#8~28 genomic mut bar
for(i in c("EGFR","MET","ERBB2","KRAS")){
  barplot.mut(gene = i)
}

for(i in setdiff(gene.highlight, c("EGFR","MET","ERBB2","KRAS"))){
  barplot.mut(gene = i, plot.mut.type = FALSE)
}

#29-31 gene fusion bar
for(i in 1:nrow(fusion.label)){
  barplot(rep(1,nrow(ITH2.table)), axes = FALSE, space = space.brp, 
          border = ifelse(ITH2.table$RNA.subtype == "No RNA sample", "gray", NA), 
          col = color.fusion[fusion.label[i,]])
  mtext(rownames(fusion.label)[i], side=2, line = 3, cex = 1, las = 2, outer = F, col = 'black')
  mtext("2.1 %", side=2, line = -2, cex = 1, las = 2, outer = FALSE, col = 'black')
}

#32 label patient ID
tmp <- barplot(rep(1,nrow(ITH2.table)), border = NA, axes = F, space = space.brp, col = "white")
text(x = sapply(unique(ITH2.table$Patient_ID), function(x){mean(tmp[ITH2.table$Patient_ID==x])}), 
     y = 0.5, label = unique(ITH2.table$Patient_ID), cex = 1, srt=90, col = "black", )
mtext("Patient ID", side=2, line = 0, cex = 1, las = 2, outer = FALSE, col = 'black')


# clonality barplot
mut.catalog <- mutate(ITH2.mut.cat, CHECK = str_glue("{Patient_ID}:{Gene}")) %>% 
  subset(!duplicated(CHECK))

tr.rate <- count(mut.catalog, Driver_gene, Clonality) %>% spread(Clonality, n)

clonality.table <- subset(mut.catalog, Gene %in% ext.data.driver.list) %>% 
  count(Gene, Clonality) %>% spread(key = Clonality, value = n, fill = 0)
clonality.fusion <- subset(ITH2.gene.fusion, Fusion_gene %in% fusion.highlight) %>% 
  mutate(CHECK = str_glue("{Patient_ID}:{Fusion_gene}"))
clonality.fusion <- clonality.fusion[duplicated(clonality.fusion$CHECK)==FALSE,] %>% 
  count(Fusion_gene, Clonality) %>% spread(key = Clonality, value = n, fill = 0)
colnames(clonality.fusion)[1] <- "Gene"

clonality.table <- rbind(clonality.table, clonality.fusion)

# perform Fisher's exact test on a contingency table against
# Trunk/Branch ratio of all non-driver genes
clonality.table$p.value <- sapply(1:nrow(clonality.table), function(x){
  cont.table <- matrix(as.numeric(c(clonality.table[x,2:3], tr.rate[2,2:3])), nrow = 2)
  p <- fisher.test(x = cont.table, alternative = "two.sided")$p.value
  return(p)
})
clonality.table$adj.pvalue <- p.adjust(clonality.table$p.value, method = "BH")

clonality.table$pvalue.mark <- case_when(
  clonality.table$adj.pvalue < 0.0001 ~ "****", 
  clonality.table$adj.pvalue < 0.001 & clonality.table$adj.pvalue >= 0.0001 ~ "***",
  clonality.table$adj.pvalue < 0.01 & clonality.table$adj.pvalue >= 0.001 ~ "**",
  clonality.table$adj.pvalue < 0.1 & clonality.table$adj.pvalue >= 0.01 ~ "*",
  TRUE ~ "")
clonality.table <- subset(clonality.table, Gene %in% c(gene.highlight, fusion.highlight)) %>% 
  arrange(desc(Truncal))
rownames(clonality.table) <- clonality.table$Gene

# gene mutation clonality barplot in transverse position
layout(matrix(1, ncol = 1))
par(mar=c(1,8,3,3))
y <- barplot(as.matrix(t(clonality.table[rev(c(gene.highlight, fusion.highlight)),3:2])), 
             col = c("blue4","steelblue1"), border = NA, las = 1, horiz = T, axes = FALSE, 
             xlim = c(0,30))
axis(side = 3, labels = seq(0,40,10), at = seq(0,40,10))
text(x = rowSums(clonality.table[rev(c(gene.highlight, fusion.highlight)),2:3]), y, pos = 4, cex = 1, las = 0, 
     label = clonality.table[rev(c(gene.highlight, fusion.highlight)), "pvalue.mark"], col = "black")
