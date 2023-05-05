### R codes for statistical testing for SBS4 activities
### across three groups and figure 3B, 3C and 3D

# Install dependency packages and load in data
rm(list = ls())
library(dplyr)
library(Rtsne)
library(ggpubr)
load("ITH2_data.Rdata")

#load in SBS activity data
ITH2.sparse <- rbind(ITH2.sparse, colSums(ITH2.sparse))
rownames(ITH2.sparse)[11] <- "Total.mut"
ITH2.sparse <- t(ITH2.sparse) %>% data.frame()
ITH2.sparse$Sector_WES_ID <- rownames(ITH2.sparse)
ITH2.sparse <- left_join(
  ITH2.sparse, dplyr::select(ITH2.table, Patient_ID, Sector_WES_ID, Group), by = "Sector_WES_ID") %>% 
  left_join(dplyr::select(ITH2.clinical, Patient_ID, Smoking_pk_yr), by = "Patient_ID")

## SBS4 assignment and z test 
SBS4.prop.by.group <- group_by(ITH2.sparse, Group) %>% 
  summarise(SBS4.assigned = sum(SBS4>0), SBS4.unassigned = sum(SBS4==0)) %>% 
  data.frame()

SBS4.prop.by.group <- cbind(SBS4.prop.by.group, 
                            SBS4.prop.by.group[,2:3]/rowSums(SBS4.prop.by.group[,2:3]))
colnames(SBS4.prop.by.group)[4:5] <- c("SBS4.assigned.prop","SBS4.unassigned.prop")


SBS4.prop.test <- lapply(list(1:2,2:3,c(1,3)), function(x){
  p  <- prop.test(as.matrix(SBS4.prop.by.group[x,2:3]))$p.value %>% as.numeric()
  df <- data.frame(Group = paste0(SBS4.prop.by.group$Group[x], collapse = " vs. "),
                   p.value = p)
  return(df)
})
SBS4.prop.test <- do.call("rbind", SBS4.prop.test)

SBS4.prop.test$adj.p.value <- p.adjust(SBS4.prop.test$p.value, method = "BH")

barplot(t(SBS4.prop.by.group[,4:5]), names.arg = SBS4.prop.by.group$driver.group, axes = FALSE)
axis(side = 2, at = seq(0,1,0.25), line = 0, las = 1, labels = paste(seq(0,100,25),"%"))

# plot the Fig 2C SBS4 activity
ITH2.SBS4 <- subset(ITH2.sparse, SBS4>0)

ggplot(ITH2.SBS4, aes(x=log10(Total.mut), y=SBS4)) + 
  geom_point(aes(color = Group), size = 3, show.legend = F) + 
  facet_grid(. ~ Group) + theme_bw() + 
  theme(axis.title.x = element_blank()) + 
  scale_color_manual(values = c("Oncogene-driven non-smoking"="#66CCEE",
                                "Oncogene-driven smoking"="#EE6677",
                                "Typical smoking"="#228833"))

# run tsne dimension reduction
tsne_out <- Rtsne(t(data.frame(ITH2.spectra.sbs96)), 
                  dims = 2, perplexity=30, verbose=TRUE, max_iter = 10000)

tmp <- data.frame(tsne_out$Y, Sector_WES_ID = colnames(ITH2.spectra.sbs96)) %>% 
  left_join(dplyr::select(ITH2.table, Patient_ID, Sector_WES_ID, Group),
            by = "Sector_WES_ID")

# plot the Fig 2D (the plot might be different from the paper)
ggplot(tmp, aes(x=X1, y=X2)) + 
  geom_point(aes(color = Group), size = 3) + 
  xlab("tSNE_1") + ylab("tSNE_2") + theme_bw() + 
  scale_color_manual(values = c("Oncogene-driven non-smoking"="#66CCEE",
                                "Oncogene-driven smoking"="#EE6677",
                                "Typical smoking"="#228833"))
