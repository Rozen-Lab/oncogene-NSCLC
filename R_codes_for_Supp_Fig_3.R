### R codes for somatic copy number alteration (SCNA)
### and Supp Figure 3A.
rm(list = ls()) 
library(tidyverse)
library(ggpubr)

load("ITH2_data.Rdata")
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
