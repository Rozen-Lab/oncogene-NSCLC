### R codes for Fig.1B & 1C (violin plots and gene enrichment plot)
# get dependency packages and load in data
rm(list = ls())
library(dplyr)
library(ggrepel)
library(ggpubr)
library(gridExtra)
load("ITH2_data.Rdata")

# compare TMB, Truncal mut, driver mut and ITH across 3 groups
ITH2.mut.cat <- mutate(ITH2.mut.cat, Gene.mut = str_glue("{Gene}_{Amino_acid_change}"))

tmp <- count(ITH2.mut.cat, Sector_ID, Clonality) %>% 
  spread(Clonality,n,fill=0) %>% 
  mutate(TMB.sector = Branch + Truncal, ITH.sector = Branch / (Branch + Truncal)) %>% 
  mutate(Branch.prop = Branch / TMB.sector, Truncal.prop = Truncal / TMB.sector) %>% 
  left_join(y = ITH2.table[,c("Patient_ID","Sector_WES_ID")], by = c("Sector_ID"="Sector_WES_ID"))

tmp2 <- group_by(tmp, Patient_ID) %>% 
  summarize(TMB.tumor = mean(TMB.sector), ITH.tumor = mean(ITH.sector)) %>% data.frame()

rownames(tmp2) <- tmp2$Patient_ID

ITH2.clinical$Tumor_mut_burden <- tmp2[ITH2.clinical$Patient_ID, "TMB.tumor"]
ITH2.clinical$ITH <- tmp2[ITH2.clinical$Patient_ID, "ITH.tumor"]

p <- lapply(c("Tumor_mut_burden","Truncal_mut","Num_driver_mut","ITH"), function(x){
  p <- ggpubr::ggviolin(
    ITH2.clinical, x = "Group", y = x, fill = "Group",
    palette =  c("#66CCEE","#EE6677","#228833"), 
    add = "boxplot", xlab = FALSE, add.params = list(fill = "white")) + 
    rremove("legend") + 
    theme(axis.title.y = element_blank(), axis.text.x = element_blank()) 
  return(p)
})
do.call("grid.arrange", c(p,ncol = 4))

# calculate odds ratio of driver mutation and p value using Fisher's test
mut.catalog <- left_join(ITH2.mut.cat, dplyr::select(ITH2.table, Sector_WES_ID, Group), 
                         by = c("Sector_ID" = "Sector_WES_ID")) %>% 
  mutate(CHECK = str_glue("{Patient_ID}:{Gene}")) %>% 
  subset(Driver_gene == "Driver" & !duplicated(CHECK)) 

gene.OR <- dplyr::count(mut.catalog, Gene, Group) %>% spread(Group, n, fill = 0)
colnames(gene.OR)[2:4] <- c("Group.1","Group.2","Group.3")
Group.n <- table(ITH2.clinical$Group)
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

# plot gene mut enrichment plot (Fig 1C)
p1 <- ggplot(gene.OR, aes(x = OR.3v1.log2, y = -log10(pvalue.3v1))) + 
  geom_point(size = 3, show.legend = FALSE,
             color = ifelse(gene.OR$pvalue.3v1 < 0.05 & gene.OR$adj.p.3v1 < 0.1, "orangered1", 
                            ifelse(gene.OR$pvalue.3v1 >= 0.05, "gray60", "steelblue3"))) + 
  geom_hline(yintercept = -log10(0.05), color = "black", linetype="dashed") + 
  ggrepel::geom_label_repel(size = 4, label.padding = 0.1, 
                            label = ifelse(gene.OR$pvalue.3v1 < 0.05, gene.OR$Gene,"")) +
  xlim(-6,6) + theme_bw() + theme(axis.title = element_blank()) + ylim(0,4.5)

p2 <- ggplot(gene.OR, aes(x = OR.3v2.log2, y = -log10(pvalue.3v2))) + 
  geom_point(size = 3, show.legend = FALSE,
             color = ifelse(gene.OR$pvalue.3v2 < 0.05 & gene.OR$adj.p.3v2 < 0.1, "orangered1", 
                            ifelse(gene.OR$pvalue.3v2 >= 0.05, "gray60", "steelblue3"))) + 
  geom_hline(yintercept = -log10(0.05), color = "black", linetype="dashed") + 
  ggrepel::geom_label_repel(size = 4, label.padding = 0.1, 
                            label = ifelse(gene.OR$pvalue.3v2 < 0.05, gene.OR$Gene,"")) +
  xlim(-6,6) + theme_bw() + theme(axis.title = element_blank()) + ylim(0,4.5) 

p3 <- ggplot(gene.OR, aes(x = OR.2v1.log2, y = -log10(pvalue.2v1))) + 
  geom_point(size = 3, show.legend = FALSE, color = "gray60") + 
  geom_hline(yintercept = -log10(0.05), color = "black", linetype="dashed") + 
  xlim(-6,6) + theme_bw() + theme(axis.title = element_blank()) + ylim(0,4.5)

ggpubr::ggarrange(p1,p2,p3,ncol=3,align="h")
