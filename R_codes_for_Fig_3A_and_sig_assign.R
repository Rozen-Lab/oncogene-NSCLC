### R codes for signature presence test, signature assignment and 
# Fig 3A (signature attribution). Courtesy of Nanhai Jiang for the 
# signature presence test and sparse assignment. 
rm(list = ls())

# Install dependency packages and load in data
library(remotes)
library(dplyr)
library(data.table)
library(parallel)

# to install mSigAct R package, please use: 
# remotes::install_github(repo = "steverozen/mSigAct", ref = "v2.3.2-branch")
library(mSigAct) 

# to install ICAMS R package, please use: 
# remotes::install_github(repo = "steverozen/ICAMS", ref = "v3.0.6-branch")
library(ICAMS) 

# to install cosmicsig R package, please use: 
# remotes::install_github(repo = "Rozen-Lab/cosmicsig", ref = "v1.0.7-branch")
library(cosmicsig) 
load("ITH2_data.Rdata")

# Get sample ids for different groups
# Group 1 are oncogene-driven non-smoking tumors, group 2 are oncogene-driven 
# smoking tumors, group 3 are typical smoking tumors
sample_ids <- list()
groups <- unique(ITH2.table$Group)
for (i in 1:3) {
  sample_ids[[i]] <- ITH2.table$Sector_WES_ID[
    ITH2.table$Group == groups[i] & ITH2.table$Tissue_type == "Tumor"]
}

sigs_sbs96_genome <- cosmicsig::COSMIC_v3.2$signature$GRCh37$SBS96

# Transform genome signatures to exome signatures
sigs_sbs96_exome <-
  ICAMS::TransformCatalog(
    catalog = sigs_sbs96_genome,
    target.region = "exome",
    target.ref.genome = "GRCh38"
  )

# Select signatures that were previously found in LUAD tumors 
# except SBS3 (HR deficiency) 
sigs_prop_sbs96 <-
  mSigAct::ExposureProportions(
    mutation.type = "SBS96",
    cancer.type = "Lung-AdenoCA"
  )

sigs_prop_sbs96 <- sigs_prop_sbs96[-3]  # exclude SBS3 from the universe

sigs_sbs96_exome_lung <-
  sigs_sbs96_exome[, names(sigs_prop_sbs96), drop = FALSE]

attribution_results <- list()
output_dir <- "./output/sig_attribution"


# Within each group, do signature presence test and sparse assignment.
# This function will probably takes 30-60 minutes if running on single core. 
for (i in 1:3) {
  sample_ids_one_group <- sample_ids[[i]]
  spectra_one_group <- ITH2.spectra.sbs96[, sample_ids_one_group, drop = FALSE]

  # Perform signature presence test to see whether SBS4 (tobacco smoking
  # signature) is existing in the tumor spectra
  tests <- mSigAct::SignaturePresenceTest(
    spectra = spectra_one_group,
    sigs = sigs_sbs96_exome_lung,
    target.sig.index = "SBS4",
    seed = 3812,
    mc.cores = 1  # can set at higher value to parellel the test
  )
  p_values <- sapply(tests, FUN = "[[", 4)

  # Do multiple testing correction
  q_values <- stats::p.adjust(p = p_values, method = "BH")

  tumor_names <- list()
  sbs4_infos <- c("no_sbs4", "sbs4")

  # Those tumors with q_value >= 0.05 do not need SBS4 in the reconstruction
  tumor_names[["no_sbs4"]] <- names(which(q_values >= 0.05))
  
  # Those tumors with q_value < 0.05 need SBS4 in the reconstruction
  tumor_names[["sbs4"]] <- names(which(q_values < 0.05))

  # Do signature attribution according to whether tumors need SBS4 in the
  # reconstruction
  for (info in sbs4_infos) {
    sample_names <- tumor_names[[info]]
    if (info == "no_sbs4") {
      sigs_to_use <-
        sigs_sbs96_exome_lung[, colnames(sigs_sbs96_exome_lung) != "SBS4",
          drop = FALSE
        ]
    } else {
      sigs_to_use <- sigs_sbs96_exome_lung
    }

    if (length(sample_names) > 0) {
      spectra_to_use <- spectra_one_group[, sample_names, drop = FALSE]
      identifier <- paste0("group_",i,"_",info)
      retval <-
        mSigAct::SparseAssignActivity(
          spectra = spectra_to_use,
          sigs = sigs_to_use,
          output.dir = file.path(output_dir, identifier),
          max.level = ncol(sigs_to_use) - 1,
          p.thresh = 0.05 / ncol(sigs_to_use),
          num.parallel.samples = 6,
          mc.cores.per.sample = 30,
          seed = 3812,
          max.subsets = .Machine$double.xmax,
          drop.low.mut.samples = FALSE
        )
      attribution_results[[identifier]] <- retval
    }
  }
}

saveRDS(object = attribution_results, 
        file = file.path(output_dir, "attribution_results.Rds"))

exposure_list <- lapply(attribution_results, FUN = function(x) {
  return(x$proposed.assignment)
})

ITH2.sparse.out <- mSigAct:::MergeListOfExposures(exposure_list)
mSigAct::WriteExposure(exposure = ITH2.sparse.out, 
                       file = file.path(output_dir, "ITH2_exposure.csv"))


#prepare clinical data to plot the signature activity plot
ITH2.sparse <- rbind(ITH2.sparse, colSums(ITH2.sparse))
rownames(ITH2.sparse)[11] <- "Total.mut"
ITH2.sparse <- t(ITH2.sparse) %>% data.frame() %>% mutate(APOBEC = SBS2 + SBS13)
ITH2.sparse$Sector_WES_ID <- rownames(ITH2.sparse)
ITH2.sparse <- left_join(ITH2.sparse, 
                         dplyr::select(ITH2.table, Group, Patient_ID, Sector_WES_ID), 
                         by = "Sector_WES_ID") %>% 
  left_join(dplyr::select(ITH2.clinical, Patient_ID, Smoking_pk_yr, Oncogene_mut), 
            by = "Patient_ID")

# prepare the color key of mut sig
color.mutsig <- c(scales::brewer_pal(palette = "Paired")(7),"grey10")
names(color.mutsig) <- c("SBS5","SBS1","SBS40","SBS18","SBS9","APOBEC","SBS17a","SBS4")

tmp <- group_by(ITH2.sparse, Patient_ID) %>% 
  summarize(mean.TMB = mean(Total.mut)) %>% data.frame()
ITH2.sparse <- left_join(ITH2.sparse, tmp, by = "Patient_ID") %>% 
  arrange(Group, desc(mean.TMB), desc(Total.mut))
ITH2.sparse.prop <- ITH2.sparse[,names(color.mutsig)]/ITH2.sparse$Total.mut

# prepare the barplot space width
space.brp <- sapply(unique(ITH2.sparse$Group), function(x){
  tmp <- subset(ITH2.sparse, Group == x)$Patient_ID %>% table()
  tmp <- tmp[unique(subset(ITH2.sparse, Group == x)$Patient_ID)]
  brp <- sapply(tmp, function(x){c(0.5,rep(0,(x-1)))}, simplify = T) %>% unlist()
  brp[1] <- 4
  return(brp)
}) %>% unlist()

color.oncogene <- c("orchid3","lightskyblue","yellowgreen","brown","tan1","gray95")
names(color.oncogene) <- c("EGFR","KRAS", "MET", "ALK", "ERBB2", "Wild type")

color.smoking <- function(x){
  case_when(x == 0 ~ "grey95",
            x > 0   & x < 20 ~ "grey75",
            x >= 20 & x < 40 ~ "grey55",
            x >= 40 & x < 60 ~ "grey35",
            x >= 60 & x < 80 ~ "grey15",
            x >= 80 ~ "grey5")
}

# plot SBS activity (Fig 2)
layout(matrix(c(rep(1,4),rep(2,10),rep(3,8),4,5,6,7,7,8,8,8), ncol = 1))
par(mar = c(0.2, 0.2, 0.1, 0.2), oma = c(1, 10, 1, 1))

#1 plot Total mut in SBS (upper)
barplot(t(ITH2.sparse[,names(color.mutsig)]), border = NA, ylim = c(2000,2500), 
        axes = FALSE, space = space.brp, axisnames = F, col = color.mutsig)
axis(side = 2, at = c(2000,2500), line = -1, las = 1, labels = c(2000,2500))

#2 plot Total mut in SBS (lower)
barplot(t(ITH2.sparse[,names(color.mutsig)]), border = NA, ylim = c(0,1250), density = NULL,
        axes = FALSE, space = space.brp, axisnames = F, col = color.mutsig)
axis(side = 2, at = c(0, 500,1000,1200), line = -1, las = 1, labels = c(0,500,1000,1200))
mtext("Signature\nactivity\n(counts)", side=2, line = 2, cex = 1, las = 2, outer = FALSE, col = 'black')

#3 plot SBS96 sig reconstruction in proportion
barplot(t(ITH2.sparse.prop[,names(color.mutsig)]), border = NA, ylim = c(0,1.05), 
        axes = FALSE, space = space.brp, axisnames = F, col = color.mutsig)
axis(side = 2, at = c(0,0.5,1), line = -1, las = 1, labels = c(0,0.5,1))
mtext("Signature\nactivity\n(proportion)", side=2, line = 2, cex = 1, 
      las = 2, outer = FALSE, col = 'black')

#4 plot smoking bar
barplot(rep(1,nrow(ITH2.sparse)), border = NA, axes = FALSE, space = space.brp, 
        col = color.smoking(ITH2.sparse$Smoking_pk_yr))
mtext("Smoking", side=2, line = -1, cex = 1, las = 2, outer = FALSE, col = 'black')

#5 plot Oncogene mut 
barplot(rep(1,nrow(ITH2.sparse)), border = NA, axes = FALSE, space = space.brp, 
        col = color.oncogene[ITH2.sparse$Oncogene_mut])
mtext("Oncogene mutation", side=2, line = -1, cex = 1, las = 2, outer = FALSE, col = 'black')

#6 plot coconut tree phylo pattern
barplot(rep(1,nrow(ITH2.sparse)), border = NA, axes = FALSE, space = space.brp, 
        col = ifelse(ITH2.sparse$Patient %in% c("A306","A590","A301","A247","A250"),
                     "grey5", "grey95"))
mtext("Phylogenic pattern", side=2, line = -1, cex = 1, las = 2, outer = FALSE, col = 'black')

#7 plot Patient ID
tmp <- barplot(rep(1,nrow(ITH2.sparse)), border = NA, axes = F, space = space.brp, col = "white")
text(x = sapply(unique(ITH2.sparse$Patient), function(x){mean(tmp[ITH2.sparse$Patient==x])}), 
     y = 0.5, label = unique(ITH2.sparse$Patient), cex = 1, srt=90, col = "black", )
mtext("Patient ID", side=2, line = -1, cex = 1, las = 2, outer = FALSE, col = 'black')

#8 plot legend
barplot(1, border = NA, axes = FALSE, col = "white", axisnames = F)
legend(x = "top", legend = names(color.mutsig), ncol=5, border = NA, bty = "n", cex = 1,
       fill = color.mutsig)
