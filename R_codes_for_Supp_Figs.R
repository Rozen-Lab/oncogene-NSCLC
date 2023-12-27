rm(list = ls()) 
library(tidyverse)
library(ggpubr)

load("ITH2_data.Rdata")

###  Supp Figure S1  ###
ITH2.cov <- read.csv("Table.S4.csv")
ITH2.clinical <- read.csv("Table.S2.csv")
ITH2.cov <- left_join(ITH2.cov, ITH2.clinical, by = "Patient_ID")
ITH2.cov$Group <- case_when(
  ITH2.cov$Group == "Non-smokers" & ITH2.cov$Mutated.oncogene == "Wild Type" ~ "NSRO-neg non-smoking", 
  ITH2.cov$Group == "Non-smokers" & ITH2.cov$Mutated.oncogene != "Wild Type" ~ "NSRO-driven non-smoking",
  ITH2.cov$Group == "NSRO-driven in smokers" ~ "NSRO-driven in smokers",
  ITH2.cov$Group == "Typical smoking" ~ "Typical smoking")
ITH2.cov$Group <- factor(ITH2.cov$Group, 
                           levels = c("NSRO-driven non-smoking", "NSRO-neg non-smoking",
                                      "NSRO-driven in smokers","Typical smoking"))
ITH2.cov$Tumor_purity <- ifelse(ITH2.cov$Tumor_purity == "< 0.1", 0, as.numeric(ITH2.cov$Tumor_purity))
ITH2.cov$Num_variant_called[is.na(ITH2.cov$Num_variant_called)] <- 0
tmp <- subset(ITH2.cov, Tissue_type == "Lung tumor") %>% 
  group_by(Patient_ID) %>% summarize(mean.var = mean(Num_variant_called))
tmp$mean.var[tmp$Patient_ID == "A037"] <- 122
ITH2.cov <- left_join(ITH2.cov, tmp, by = "Patient_ID") %>% 
  arrange(Group, desc(mean.var), desc(Tissue_type), desc(Num_variant_called))

space.brp <- sapply(table(ITH2.cov$Patient_ID)[unique(ITH2.cov$Patient_ID)], function(x){
  c(1,rep(0,(x-1)))}, simplify = T) %>% unlist()

color.driver <- c("EGFR"="orchid3","KRAS"="lightskyblue","BRAF"="steelblue4",
                  "MET"="yellowgreen","HER2"="tan1","ALK"="brown","Wild Type"="gray90")

layout(matrix(c(1,1,1,1,2,2,3,3,4,5,6), ncol = 1))
par(mar = c(0.2, 0.2, 0.1, 0.2), oma = c(1, 2, 2, 0.5))

barplot(ITH2.cov$Num_variant_called, axes = FALSE, space = space.brp, 
        col = ifelse(ITH2.cov$Tumor_purity == 0, "red", "gray10"), border = NA)
axis(side = 2, at = seq(0,2500,500), labels = seq(0,2500,500), line = -1, las = 1)
mtext("Num\nVar", side=2, line = 4, cex = 1, las = 2, outer = F, col = 'black')

barplot(matrix(c(ITH2.cov$Tumor_purity, 1-ITH2.cov$Tumor_purity), nrow=2, byrow=T), border = NA, 
        axes = FALSE, col = c("gray10","gray80"), space = space.brp, ylim = c(0,1.1))
abline(h=0.1, lty=3, lwd=1, col="orangered")
axis(side = 2, at = c(0,0.5,1), labels = c(0,0.5,1), line = -1, las = 1)
mtext("Tumor\npurity", side=2, line = 4, cex = 1, las = 2, outer = F, col = 'black')

barplot(ITH2.cov$Mean_cov_per_target, border = NA, axes = FALSE, ylim=c(0,110),
        col = "gray10", space = space.brp)
axis(side = 2, at = seq(0,100,25), labels = seq(0,100,25), line = -1, las = 1)

barplot(rep(1,nrow(ITH2.cov)), border = NA, axes = FALSE, space = space.brp, 
        col = case_when(
          ITH2.cov$Tissue_type == "Lung tumor" & ITH2.cov$Tumor_purity >= 0.1 ~ "steelblue1",
          ITH2.cov$Tissue_type == "Lung tumor" & ITH2.cov$Tumor_purity == 0 ~ "gray",
          ITH2.cov$Tissue_type != "Lung tumor" ~ "plum1",
        ))
mtext("Tissue type", side=2, line = -1, cex = 1, las = 2, outer = FALSE, col = 'black')

barplot(rep(1,nrow(ITH2.cov)), border = NA, axes = FALSE, space = space.brp, 
        col = color.driver[ITH2.cov$Mutated.oncogene])
mtext("Oncogene mutation", side=2, line = -1, cex = 1, las = 2, outer = FALSE, col = 'black')

tmp <- barplot(rep(1,nrow(ITH2.cov)), border = NA, axes = F, space = space.brp, col = "white")
text(x = sapply(unique(ITH2.cov$Patient_ID), function(x){mean(tmp[ITH2.cov$Patient_ID==x])}), 
     y = 0.5, label = unique(ITH2.cov$Patient_ID), cex = 1, srt=90, col = "black", )
mtext("Patient ID", side=2, line = -1, cex = 1, las = 2, outer = FALSE, col = 'black')



### Supp Figure S2 ###
GIS.maf <- read.table("data-GIS_snv_indel.maf", sep = "\t", header = T) %>% 
  subset(Chromosome %in% c(1:22, "X","Y") & Tumor_Sample_Barcode %in% GIS.clinical$Patient.ID &
           Variant_Classification %in% 
           c("Missense_Mutation","Nonsense_Mutation","Splice_Site","In_Frame_Del",
             "In_Frame_Ins","Frame_Shift_Del", "Frame_Shift_Ins"))

GIS.clinical <- read.csv("data-GIS_clinical.csv")
GIS.clinical$driver.type <- sapply(1:nrow(GIS.clinical), function(x){
  if(GIS.clinical$driver.gene[x] %in% c("Wild Type", "ALK")){
    return(GIS.clinical$driver.gene[x])
  } else if(GIS.clinical$driver.gene[x] == "EGFR"){
    maf <- subset(GIS.maf.NS, Hugo_Symbol == "EGFR" & Tumor_Sample_Barcode == GIS.clinical$Patient.ID[x] &
                    !HGVSp_Short %in% c("p.L62R","p.D256Y","p.R973L","p.A566S","p.E548Q")) %>% 
      mutate(AAChange = sub("p.","",HGVSp_Short))
    return(paste(maf$AAChange, collapse = ";"))
  } else if(GIS.clinical$driver.gene[x] == "ERBB2"){
    maf <- subset(GIS.maf.NS, Hugo_Symbol == "ERBB2" & Tumor_Sample_Barcode == GIS.clinical$Patient.ID[x] &
                    Variant_Classification == "In_Frame_Ins") %>% 
      mutate(AAChange = sub("p.","",HGVSp_Short))
    return(maf$AAChange)
  } else if(GIS.clinical$driver.gene[x] == "MET"){
    maf <- subset(GIS.maf, Hugo_Symbol == "MET" & Tumor_Sample_Barcode == GIS.clinical$Patient.ID[x] & 
                    Variant_Classification %in% c("Splice_Site","Intron"))
    return(maf$HGVSc)
  } else if(GIS.clinical$driver.gene[x] %in% c("KRAS","BRAF","KRAS.BRAF")){
    maf <- subset(GIS.maf.NS, Hugo_Symbol %in% c("KRAS","BRAF") & 
                    Tumor_Sample_Barcode == GIS.clinical$Patient.ID[x]) %>% 
      mutate(AAChange = sub("p.","",HGVSp_Short))
    return(paste(maf$AAChange, collapse = ";"))
  } 
})

GIS.clinical <- subset(GIS.clinical, driver.group != "Group 4")
driver.list <- readRDS("data-driver.list.COSMIC.ext.rds")
GIS.maf$driver <- ifelse(GIS.maf$Hugo_Symbol %in% driver.list, "Driver", "Non-driver")
tmp <- subset(GIS.maf, driver == "Driver") %>% count(Tumor_Sample_Barcode)
colnames(tmp)[2] <- "Driver.n"
GIS.clinical <- left_join(GIS.clinical, tmp, by = c("Patient.ID"="Tumor_Sample_Barcode"))
GIS.clinical$Driver.n[is.na(GIS.clinical$Driver.n)] <- 0

compare.group <- list(c("Group 1","Group 2"), c("Group 2","Group 3"), c("Group 1","Group 3"))
p <- lapply(c("Mutation.Count","Driver.n"), function(x){
  y.max <- case_when(
    x == "Mutation.Count" ~ 1500,
    x == "Driver.n" ~ 120
  )
  p <- ggviolin(GIS.clinical, x = "driver.group", 
                y = x, fill = "driver.group", 
                palette =  c("#66CCEE","#EE6677","#228833"), width=0.8, 
                add = "boxplot", xlab = FALSE, add.params = list(fill = "white")) + rremove("legend") + 
    theme(axis.title.y = element_blank(), axis.text.x = element_blank()) + 
    scale_y_continuous(limits = c(0,y.max)) + 
    stat_compare_means(comparisons = compare.group, label = "p.format", hide.ns = F, show.legend = F)
  return(p)
})
do.call("ggarrange", c(p,ncol = 2))

GIS.driver <- subset(GIS.maf, driver == "Driver") %>% 
  left_join(y = dplyr::select(GIS.clinical, Patient.ID, driver.group), 
            by = c("Tumor_Sample_Barcode" = "Patient.ID")) %>% 
  mutate(CHECK = str_glue("{Tumor_Sample_Barcode}:{Hugo_Symbol}")) %>% 
  subset(!duplicated(CHECK))

gene.OR <- dplyr::count(GIS.driver, Hugo_Symbol, driver.group) %>% 
  spread(driver.group, n, fill = 0)
colnames(gene.OR)[2:4] <- c("Group.1","Group.2","Group.3")
Group.n <- table(GIS.clinical$driver.group)
gene.OR <- mutate(gene.OR, Group.1.wt = Group.n[1]-Group.1, 
                  Group.2.wt = Group.n[2]-Group.2, Group.3.wt = Group.n[3]-Group.3)
gene.OR$prev <- rowSums(gene.OR[2:4])*100/sum(Group.n)

tmp <- lapply(1:nrow(gene.OR), function(x){
  tab.OR <- data.frame(MT = as.numeric(gene.OR[x,2:4]), 
                       MT.1 = as.numeric(gene.OR[x,2:4])+1,
                       WT = as.numeric(gene.OR[x,5:7]))
  OR.2v1 <- (tab.OR[2,1]/tab.OR[2,3])/(tab.OR[1,1]/tab.OR[1,3])
  OR.2v1.plus1 <- (tab.OR[2,2]/tab.OR[2,3])/(tab.OR[1,2]/tab.OR[1,3])
  pvalue.2v1 <- fisher.test(as.matrix(tab.OR[c(1,2),c(1,3)]))$p.value
  OR.3v1 <- (tab.OR[3,1]/tab.OR[3,3])/(tab.OR[1,1]/tab.OR[1,3])
  OR.3v1.plus1 <- (tab.OR[3,2]/tab.OR[3,3])/(tab.OR[1,2]/tab.OR[1,3])
  pvalue.3v1 <- fisher.test(as.matrix(tab.OR[c(1,3),c(1,3)]))$p.value
  OR.3v2 <- (tab.OR[3,1]/tab.OR[3,3])/(tab.OR[2,1]/tab.OR[2,3])
  OR.3v2.plus1 <- (tab.OR[3,2]/tab.OR[3,3])/(tab.OR[2,2]/tab.OR[2,3])
  pvalue.3v2 <- fisher.test(as.matrix(tab.OR[c(2,3),c(1,3)]))$p.value
  
  return(data.frame(OR.2v1, OR.2v1.plus1, pvalue.2v1, 
                    OR.3v1, OR.3v1.plus1, pvalue.3v1, 
                    OR.3v2, OR.3v2.plus1, pvalue.3v2))
})
tmp <- do.call("rbind",tmp)
gene.OR <- cbind(gene.OR, tmp)
gene.OR <- mutate(gene.OR, 
                  OR.2v1.log2 = log2(OR.2v1.plus1), 
                  adj.p.2v1 = p.adjust(pvalue.2v1, method = "BH"),
                  OR.3v1.log2 = log2(OR.3v1.plus1), 
                  adj.p.3v1 = p.adjust(pvalue.3v1, method = "BH"),
                  OR.3v2.log2 = log2(OR.3v2.plus1), 
                  adj.p.3v2 = p.adjust(pvalue.3v2, method = "BH")
)

p1 <- ggplot(gene.OR, aes(x = OR.3v1.log2, y = -log10(pvalue.3v1))) + 
  geom_point(size = 3, show.legend = FALSE,
             color = ifelse(gene.OR$pvalue.3v1 < 0.05 & gene.OR$adj.p.3v1 < 0.1, "orangered1", 
                            ifelse(gene.OR$pvalue.3v1 >= 0.05, "gray60", "steelblue3"))) + 
  geom_hline(yintercept = -log10(0.05), color = "black", linetype="dashed") + 
  ggrepel::geom_label_repel(size = 4, label.padding = 0.1, 
                            label = ifelse(gene.OR$pvalue.3v1 < 0.05, gene.OR$Hugo_Symbol,"")) +
  xlim(-8,8) + theme_bw() + theme(axis.title = element_blank()) + ylim(0,24)

p2 <- ggplot(gene.OR, aes(x = OR.3v2.log2, y = -log10(pvalue.3v2))) + 
  geom_point(size = 3, show.legend = FALSE,
             color = ifelse(gene.OR$pvalue.3v2 < 0.05 & gene.OR$adj.p.3v2 < 0.1, "orangered1", 
                            ifelse(gene.OR$pvalue.3v2 >= 0.05, "gray60", "steelblue3"))) + 
  geom_hline(yintercept = -log10(0.05), color = "black", linetype="dashed") + 
  ggrepel::geom_label_repel(size = 4, label.padding = 0.1, 
                            label = ifelse(gene.OR$pvalue.3v2 < 0.05, gene.OR$Hugo_Symbol,"")) +
  xlim(-8,8) + theme_bw() + theme(axis.title = element_blank()) + ylim(0,24) 

p3 <- ggplot(gene.OR, aes(x = OR.2v1.log2, y = -log10(pvalue.2v1))) + 
  geom_point(size = 3, show.legend = FALSE,
             color = ifelse(gene.OR$pvalue.2v1 < 0.05 & gene.OR$adj.p.2v1 < 0.1, "orangered1", 
                            ifelse(gene.OR$pvalue.2v1 >= 0.05, "gray60", "steelblue3"))) + 
  geom_hline(yintercept = -log10(0.05), color = "black", linetype="dashed") + 
  ggrepel::geom_label_repel(size = 4, label.padding = 0.1, 
                            label = ifelse(gene.OR$pvalue.2v1 < 0.05, gene.OR$Hugo_Symbol,"")) +
  xlim(-8,8) + theme_bw() + theme(axis.title = element_blank()) + ylim(0,24)

ggpubr::ggarrange(p1,p2,p3,ncol=3,align="h")




###  Supp Figure S5  ###
ITH2.hatchet <- readcsv("data-SCNA-segment-HATCHet.csv")
ctm <- read.delim("data-GRCh38p7-region.tsv")

### binning the genome at size of 5MB and calculate mean CN
CNV.summary <- function(x, res.mb = 5){
  res <- res.mb*1000000
  cat(paste("Analyzing CNV segments of", x, "\n"))
  seg.ID <- subset(ITH2.hatchet, ID == x)

  seg.nc <- lapply(paste0("chr",1:22), function(x){
    seg.chr <- subset(seg.ID, CHR == x)
    df1 <- subset(seg.chr, END <= ctm$centromere_start[ctm$chr==x])
    df2 <- subset(seg.chr, START > ctm$centromere_end[ctm$chr==x])
    
    df3 <- subset(seg.chr, START <= ctm$centromere_start[ctm$chr==x] & 
                    END > ctm$centromere_start[ctm$chr==x] & 
                    END <= ctm$centromere_end[ctm$chr==x])
    if(nrow(df3)>0){
      df3$END <- ctm$centromere_start[ctm$chr==x]
    }
    
    df4 <- subset(seg.chr, END > ctm$centromere_end[ctm$chr==x] & 
                    START > ctm$centromere_start[ctm$chr==x] & 
                    START <= ctm$centromere_end[ctm$chr==x])
    if(nrow(df4)>0){
      df4$START <- ctm$centromere_end[ctm$chr==x]
    }
    
    df5 <- subset(seg.chr, START <= ctm$centromere_start[ctm$chr==x] & 
                    END > ctm$centromere_end[ctm$chr==x])
    if(nrow(df5)>0){
      df5 <- df5[c(1,1),]
      df5$END[1] <- ctm$centromere_start[ctm$chr==x]
      df5$START[2] <- ctm$centromere_end[ctm$chr==x]
    }
    
    df <- rbind(df1,df2,df3,df4,df5) %>% arrange(START)
    return(df)
  })
  seg.nc <- do.call("rbind", seg.nc)
  
  seg.cnv <- lapply(paste0("chr",1:22), function(x){
    seg.chr <- subset(seg.nc, CHR == x)
    max.bin <- ceiling(max(seg.chr$END)/res)
    seg.bin <- seq(0, res*max.bin, res)

    seg.cn <- lapply(1:max.bin, function(x){
      df1 <- subset(seg.chr, START > seg.bin[x] & END <= seg.bin[x+1])
      
      df2 <- subset(seg.chr, START <= seg.bin[x] & END > seg.bin[x+1])
      if(nrow(df2)>0){
        df2$START <- seg.bin[x]
        df2$END <- seg.bin[x+1]
      }
      
      df3 <- subset(seg.chr, START <= seg.bin[x] & END > seg.bin[x] & END <= seg.bin[x+1])
      if(nrow(df3)>0){
        df3$START <- seg.bin[x]
      }
      
      df4 <- subset(seg.chr, START > seg.bin[x] & START <= seg.bin[x+1] & END > seg.bin[x+1])
      if(nrow(df4)>0){
        df4$END <- seg.bin[x+1]
      }
      
      df <- rbind(df1, df2, df3, df4)
      if(nrow(df)>0){
        CNa.sum <- (sum(df$CNa*(df$END-df$START))/sum(df$END-df$START)) %>% 
          round()
        CNb.sum <- (sum(df$CNb*(df$END-df$START))/sum(df$END-df$START)) %>% 
          round()
        CNt.sum <- CNa.sum + CNb.sum
      } else{
        CNa.sum <- NA
        CNb.sum <- NA
        CNt.sum <- NA
      }
      tmp <- data.frame(START = seg.bin[x]+1, END = seg.bin[x+1], 
                        CNa.sum = CNa.sum, CNb.sum = CNb.sum, CNt.sum = CNt.sum)
      return(tmp)
    })
    seg.cn <- do.call("rbind", seg.cn)
    seg.cn$CHR <- x
    return(seg.cn)
  })
  seg.cnv <- do.call("rbind",seg.cnv)
  seg.cnv$ID <- x
  seg.cnv$BG.ploidy <- seg.ID$BG.ploidy[1]
  seg.cnv$Group <- seg.ID$Group[1]
  
  return(seg.cnv)
}

ITH2.seg <- lapply(unique(ITH2.hatchet$ID), CNV.summary)
ITH2.seg <- do.call("rbind", ITH2.seg)

ITH2.seg <- mutate(ITH2.seg, CNV = CNt.sum - BG.ploidy, 
                   START = as.integer(START), END = as.integer(END)) %>% 
  mutate(SEG = str_glue("{CHR}:{START}-{END}")) %>% 
  dplyr::select(ID, CHR, START, END, SEG, starts_with("CN"), BG.ploidy, Group)
ITH2.seg$CNa.sum[is.nan(ITH2.seg$CNa.sum)] <- NA
ITH2.seg$CNb.sum[is.nan(ITH2.seg$CNb.sum)] <- NA
ITH2.seg$CNt.sum[is.nan(ITH2.seg$CNt.sum)] <- NA
ITH2.seg$CNV[is.nan(ITH2.seg$CNV)] <- NA

ITH2.seg <- mutate(ITH2.seg, CNV.status = case_when(
  !is.na(CNt.sum) & CNa.sum == 0 & CNb.sum == 0 ~ "CN.del",
  !is.na(CNt.sum) & CNa.sum > 0 & CNV == 1 & CNb.sum > 0 ~ "CN1.gain.het",
  !is.na(CNt.sum) & CNa.sum > 0 & CNV == 1 & CNb.sum == 0 ~ "CN1.gain.loh",
  !is.na(CNt.sum) & CNa.sum > 0 & CNV >= 2 & CNb.sum > 0 ~ "CN2.gain.het",
  !is.na(CNt.sum) & CNa.sum > 0 & CNV >= 2 & CNb.sum == 0 ~ "CN2.gain.loh",
  !is.na(CNt.sum) & CNa.sum > 0 & CNV == -1 & CNb.sum > 0 ~ "CN1.loss.het",
  !is.na(CNt.sum) & CNa.sum > 0 & CNV == -1 & CNb.sum == 0 ~ "CN1.loss.loh",
  !is.na(CNt.sum) & CNa.sum > 0 & CNV <= -2 & CNb.sum > 0 ~ "CN2.loss.het",
  !is.na(CNt.sum) & CNa.sum > 0 & CNV <= -2 & CNb.sum == 0 ~ "CN2.loss.loh",
  !is.na(CNt.sum) & CNa.sum > 0 & CNV == 0 & CNb.sum == 0 ~ "CN0.loh",
  !is.na(CNt.sum) & CNa.sum > 0 & CNV == 0 & CNb.sum > 0 ~ "CN.neutral",
  TRUE ~ "Unknown"
))

CNV.sum <- data.frame(SEG = unique(ITH2.seg$SEG)) %>% 
  separate(col = "SEG", into = c("CHR","RANGE"), sep = ":", remove = F) %>% 
  separate(col = "RANGE", into = c("START","END"), sep = "-", remove = F, convert = T) %>% 
  mutate(CHR = factor(CHR, levels = paste0("chr",1:22))) %>% 
  arrange(CHR, START)
CNV.sum$pos <- 1:nrow(CNV.sum)
plot.chr <- group_by(CNV.sum, CHR) %>% 
  summarize(width = length(CHR), break.start = min(pos)-1, break.end = max(pos)) %>% 
  mutate(pos = break.end-0.5*width)

CNV.gain <- count(ITH2.seg, SEG, CNV.status) %>% 
  subset(str_detect(CNV.status, "gain")) %>% 
  mutate(CNV.status = factor(CNV.status, levels = c("CN2.gain.loh","CN2.gain.het",
                                                    "CN1.gain.loh","CN1.gain.het")),
         SEG = factor(SEG, levels = unique(CNV.sum$SEG))) %>% 
  left_join(y = CNV.sum[,c("SEG","pos")], by = "SEG")

CNV.loss <- count(ITH2.seg, SEG, CNV.status) %>% 
  subset(str_detect(CNV.status, "loss|CN0|del")) %>% 
  mutate(CNV.status = factor(CNV.status, levels = c("CN.del","CN2.loss.loh","CN2.loss.het",
                                                    "CN1.loss.loh","CN1.loss.het","CN0.loh")),
         SEG = factor(SEG, levels = unique(CNV.sum$SEG))) %>% 
  left_join(y = CNV.sum[,c("SEG","pos")], by = "SEG")

pt.n <- n_distinct(ITH2.seg$ID)

CNV.gene <- read.delim("data-CN-driver-gene.tsv")
CNV.gene$plot.pos <- sapply(1:nrow(CNV.gene), function(x){
  round(CNV.gene$Pos[x]/5000000, digits = 2) + plot.chr$break.start[plot.chr$CHR == CNV.gene$Chr[x]]
})

color.CNV <- c("CN2.gain.loh"="red3","CN2.gain.het"="salmon2",
               "CN1.gain.loh"="sienna","CN1.gain.het"="tan1",
               "CN2.loss.loh"="steelblue4","CN2.loss.het"="steelblue1",
               "CN1.loss.loh"="limegreen","CN1.loss.het"="palegreen",
               "CN0.loh"="gray50", "CN.del"="darkviolet")

p1 <- ggplot(CNV.gain) + 
  geom_rect(data=plot.chr, mapping=aes(xmin=break.start, xmax=break.end, ymin=0, ymax=80),
            fill = rep(c("gray90","white"),11)) + 
  geom_segment(data=subset(CNV.gene, CNV.type =="CN.gain"), 
               mapping=aes(x=plot.pos, xend=plot.pos, y=0, yend=70), color="black", size=0.25) + 
  ggrepel::geom_text_repel(data=subset(CNV.gene, CNV.type =="CN.gain"), 
            mapping=aes(x=plot.pos, label=Gene, y=70), size = 3) + 
  geom_hline(yintercept = seq(0,60,20), color = "gray80", size=0.25) + 
  geom_col(aes(x=pos, y=100*n/pt.n, fill=CNV.status), width=1, show.legend = F) + 
  scale_x_continuous(breaks = c(0,plot.chr$break.end), minor_breaks = NULL, labels = NULL) +
  theme_void() + theme(axis.title = element_blank()) + 
  scale_fill_manual(values = color.CNV)

p2 <- ggplot(plot.chr) + 
  geom_rect(aes(xmin=break.start, xmax=break.end, ymin=0, ymax=1), 
            fill = rep(c("gray90","white"),11), color="black") + 
  geom_text(aes(x=pos, y=0.5, label=1:22), size = 4) + 
  theme_void()

p3 <- ggplot(CNV.loss) + 
  geom_rect(data=plot.chr, mapping=aes(xmin=break.start, xmax=break.end, ymin=0, ymax=-100),
            fill = rep(c("gray90","white"),11)) + 
  geom_segment(data=subset(CNV.gene, CNV.type =="CN.loss"), 
               mapping=aes(x=plot.pos, xend=plot.pos, y=0, yend=-85), color="black", size=0.25) + 
  ggrepel::geom_text_repel(data=subset(CNV.gene, CNV.type =="CN.loss"), 
                           mapping=aes(x=plot.pos, label=Gene, y=-85), size = 3) + 
  geom_hline(yintercept = -seq(0,80,20), color = "gray80", size=0.25) + 
  geom_col(aes(x=pos, y=-100*n/pt.n, fill=CNV.status), width=1, show.legend = F) + 
  geom_segment(x=0, xend=20, y=-70, yend=-70, color="black") + 
  annotate("text",x=10, y=-72, label="100MB", size = 4) + 
  scale_x_continuous(breaks = c(0,plot.chr$break.end), minor_breaks = NULL, labels = NULL) +
  theme_void() + theme(axis.title = element_blank()) + 
  scale_fill_manual(values = color.CNV)

ggpubr::ggarrange(p1,p2,p3,heights = c(8,1,10), ncol=1, align="v")


## comparison of TMB/ITH between male vs. female NSRO-driven tumors in GIS cohort
GIS.clinical <- read.csv("ITH2/GIS cohort/GIS031_clinical_valid.csv") %>% 
  subset(driver.group %in% c("Group 1", "Group 2"))

GIS.maf <- read.table("ITH2/GIS cohort/GIS_snv_indel.maf", sep = "\t", header = T) %>% 
  subset(Chromosome %in% c(1:22, "X","Y") & Tumor_Sample_Barcode %in% GIS.clinical$Patient.ID &
           Variant_Classification %in% 
           c("Missense_Mutation","Nonsense_Mutation","Splice_Site","In_Frame_Del",
             "In_Frame_Ins","Frame_Shift_Del", "Frame_Shift_Ins"))
driver.list <- readRDS("ITH2/Driver_gene_list/driver.list.COSMIC.ext.rds")
GIS.maf$driver <- ifelse(GIS.maf$Hugo_Symbol %in% driver.list, "Driver", "Non-driver")
tmp <- subset(GIS.maf, driver == "Driver") %>% count(Tumor_Sample_Barcode)
colnames(tmp)[2] <- "Driver.n"
GIS.clinical <- left_join(GIS.clinical, tmp, by = c("Patient.ID"="Tumor_Sample_Barcode"))
GIS.clinical$Driver.n[is.na(GIS.clinical$Driver.n)] <- 0

library(ggpubr)
compare.group <- list(c("Male","Female"))
p <- lapply(c("Mutation.Count","Driver.n"), function(x){
  y.max <- case_when(
    x == "Mutation.Count" ~ 300,
    x == "Driver.n" ~ 50
  )
  p <- ggviolin(GIS.clinical, x = "Gender", 
                y = x, fill = "Gender", 
                palette = c("orchid1","turquoise4"), width=0.8, 
                add = "boxplot", xlab = FALSE, add.params = list(fill = "white")) + rremove("legend") + 
    theme(axis.title.y = element_blank(), axis.text.x = element_blank()) + 
    scale_y_continuous(limits = c(0,y.max)) + 
    stat_compare_means(comparisons = compare.group, label = "p.format", hide.ns = F, show.legend = F)
  return(p)
})
do.call("ggarrange", c(p,ncol = 2))
table(GIS.clinical$Gender)

