### R codes for classifying TRU vs. non-TRU gene expression
### subtypes, pathway activity transforming, and Fig. 4 plot

# Install dependency packages and load in data
rm(list = ls()) 
library(dplyr)
library(stringr)
library(ConsensusClusterPlus)
library(ComplexHeatmap)
library(GSVA)
library(GSEABase)
library(org.Hs.eg.db)
library(clusterProfiler)
library(limma)
load("ITH2_data.Rdata")

# classification of gene expression subtypes (TRU vs. non-TRU) using 
# predictor gene panel and methods provided by Wilkerson 2012.
# The gene panel has been loaded in as 'ext.data.TRU.gene' object.
View(ext.data.TRU.gene)
rownames(ext.data.TRU.gene) <- ext.data.TRU.gene$gene.correct

RNA.matrix <- ITH2.RNA.TPM[ext.data.TRU.gene$gene.correct,] # subset RNA matrix by the gene panel
RNA.matrix <- RNA.matrix/rowMeans(RNA.matrix)
RNA.matrix <- sweep(RNA.matrix, 1, apply(RNA.matrix,1,median,na.rm=T))

cluster.results <- ConsensusClusterPlus(
  d = as.matrix(RNA.matrix), 
  maxK = 5,
  reps = 10000,
  pItem = 0.8,
  pFeature = 1,
  title = "./TRU.cluster",
  clusterAlg = "km", 
  distance = "euclidean", 
  seed = 12620,
  plot = "pdf"
)

cluster.results[[2]]   ## we will use clustering results with k of 2 

cluster.TRU <- colnames(RNA.matrix)[cluster.results[[2]]$consensusClass == 1] ## TRU

cluster.non.TRU <- colnames(RNA.matrix)[cluster.results[[2]]$consensusClass == 2] ## Non-TRU

## visualize heatmap of TRU vs. non-TRU
Heatmap(as.matrix(RNA.matrix[,c(cluster.TRU, cluster.non.TRU)]), 
        column_split = rep(c("TRU", "Non-TRU"), 
                           c(length(cluster.TRU), length(cluster.non.TRU))), 
        column_names_gp = gpar(fontsize = 6), 
        row_names_gp = gpar(fontsize = 6), 
        cluster_rows = F) 


### transform transcript level gene expression to pathway level
### gene expression using GSVA package

ext.data.reactome <- GSEABase::getGmt(
  con = "./External data/c2.cp.reactome.v2022.1.Hs.symbols.gmt", 
  geneIdType = SymbolIdentifier(), 
  collectionType = BroadCollection(category="h"), 
  sep = "\t"
)

ITH2.RNA.gsva <- GSVA::gsva(
  expr = as.matrix(ITH2.RNA.TPM), 
  gset.idx.list = ext.data.reactome, 
  min.sz=10, max.sz=500
)


### conduct differential pathway analysis and clustering
ITH2.RNA$Smoking_status <- ifelse(ITH2.RNA$Smoking_pk_yr > 0, "Ex.Current smoker", "Never smoker")
ITH2.RNA$SBS4 <- as.numeric(ITH2.sparse["SBS4",ITH2.RNA$Sector_WES_ID])
ITH2.RNA$SBS4.activity <- ifelse(ITH2.RNA$SBS4 > median(ITH2.RNA$SBS4[ITH2.RNA$SBS4>0]), 
                                 "SBS4.high", "SBS4.low")

ITH2.RNA$Group <- factor(ITH2.RNA$Group, 
                         levels = c("Non-smoking","NSRO-driven smoking","Typical-smoking"))

mod <- model.matrix(~ ITH2.RNA$Smoking_status, levels = c("Ex.Current smoker","Never smoker"))
colnames(mod)[2] <- c("non.smoker")
fit <- lmFit(object = ITH2.RNA.gsva, design = mod)
fit <- eBayes(fit)
tt <- topTable(fit, coef=2, n=Inf, adjust.method = "fdr")
tt <- rbind(arrange(subset(tt, adj.P.Val <= 0.05 & logFC > 0), adj.P.Val), 
            arrange(subset(tt, adj.P.Val <= 0.05 & logFC < 0), adj.P.Val))

DEpwys <- rownames(tt) # differential pathways used for clustering

gsva.pathway <- ITH2.RNA.gsva[DEpwys,ITH2.RNA$Sector_RNA_ID] %>% t() %>% scale() %>% t()

hm <- ComplexHeatmap::Heatmap(as.matrix(gsva.pathway), cluster_columns = T, column_km = 2) 

ITH2.RNA$pathway.group <- "Group I" 
ITH2.RNA$pathway.group[column_order(hm)[[2]]] <- "Group II" 
ITH2.RNA$gene.exp.subtype <- ifelse(ITH2.RNA$Sector_RNA_ID %in% cluster.TRU, "TRU", "Non-TRU")

tmp <- dplyr::count(ITH2.RNA, Patient_ID, pathway.group) %>% 
  tidyr::spread(pathway.group, n, fill=0) %>% 
  mutate(pathway.group.tumor = ifelse(`Group I`==0 , "Group II", 
                                      ifelse(`Group II`==0, "Group I", "Mixed"))) %>% 
  mutate(pathway.group.tumor = factor(pathway.group.tumor, 
                                      levels = c("Group I","Mixed","Group II"))) %>% 
  left_join(ITH2.RNA[!duplicated(ITH2.RNA$Patient_ID),c("Group", "Patient_ID")], 
            by = "Patient_ID") %>% 
  arrange(Group, pathway.group.tumor, desc(`Group II`), `Group I`)

ITH2.RNA <- mutate(ITH2.RNA, Patient_ID = factor(Patient_ID, levels = tmp$Patient_ID)) %>% 
  arrange(Patient_ID, pathway.group) %>% 
  mutate(Patient.mark = rep(rep(c("A","B"),16), table(Patient_ID)))

ITH2.RNA$Smoking_pk_yr <- case_when(
  ITH2.RNA$Smoking_pk_yr == 0 ~ "0",
  ITH2.RNA$Smoking_pk_yr > 0   & ITH2.RNA$Smoking_pk_yr < 20 ~ "0-20",
  ITH2.RNA$Smoking_pk_yr >= 20 & ITH2.RNA$Smoking_pk_yr < 40 ~ "20-40",
  ITH2.RNA$Smoking_pk_yr >= 40 & ITH2.RNA$Smoking_pk_yr < 60 ~ "40-60",
  ITH2.RNA$Smoking_pk_yr >= 60 & ITH2.RNA$Smoking_pk_yr < 80 ~ "60-80",
  ITH2.RNA$Smoking_pk_yr >= 80 ~ "80up")

### plot activities of top 10 up- and down-regulated differential pathways on a heatmap 
DEpwys.up <- subset(tt, adj.P.Val <= 0.05 & logFC > 0) %>% top_n(n = 10, wt = -(adj.P.Val)) %>% 
  rownames() # differential pathways upregulated in non-smoking NSCLC

DEpwys.down <- subset(tt, adj.P.Val <= 0.05 & logFC < 0) %>% top_n(n = 10, wt = -(adj.P.Val)) %>% 
  rownames() # differential pathways upregulated in non-smoking NSCLC

pathways <- c(DEpwys.up, DEpwys.down)

gsva.pathway <- ITH2.RNA.gsva[pathways, ITH2.RNA$Sector_RNA_ID] %>% t() %>% scale() %>% t()

pathway.q.value <- tt[pathways, "adj.P.Val"] %>% format(digits = 2)

rownames(gsva.pathway) <- paste(str_pad(pathway.q.value, width = 8, side = "right"), 
                                rownames(gsva.pathway), sep = " ")

Heatmap(as.matrix(gsva.pathway), name = "Z score", 
        column_split = rep(unique(ITH2.RNA$Group), table(ITH2.RNA$Group)),
        cluster_columns = F, column_names_gp = gpar(fontsize = 4), 
        row_split = rep(c("Up-regulated","Down-regulated"), 
                        c(length(DEpwys.up), length(DEpwys.down))), 
        cluster_rows = F, row_km = 1, row_names_gp = gpar(fontsize = 6), 
        top_annotation = HeatmapAnnotation(
          Group = ITH2.RNA$Group, 
          Oncogene = ITH2.RNA$Oncogene_mut, 
          TRU.subtype = ITH2.RNA$gene.exp.subtype,
          SBS4.activity = ITH2.RNA$SBS4.activity, 
          Smoking = ITH2.RNA$Smoking_pk_yr,
          Patient.mark = ITH2.RNA$Patient.mark,
          col = list(Group = c("Non-smoking"="#66CCEE", 
                               "NSRO-driven smoking"="#EE6677",
                               "Typical-smoking"="#228833"),
                     Oncogene = c("EGFR"="orchid3","KRAS"="lightskyblue","MET"="yellowgreen",
                                  "ALK"="brown","ERBB2"="tan1","Wild type"="gray90"),
                     TRU.subtype = c("TRU"="violetred1", "Non-TRU"="steelblue1"),
                     SBS4.activity = c("SBS4.high"="black", "SBS4.low"="gray90"),
                     Patient.mark = c("A" = "gray90", "B" = "gray10"),
                     Smoking = c("0"="grey90", "0-20"="grey75","20-40"="grey55",
                                 "40-60"="grey35","60-80"="grey15","80up"="grey5")
          )
        ))
