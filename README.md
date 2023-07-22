# About this code repository
This repository provides R codes and data used for the paper

> Huang, C.-Y., et al. (2023). "Dearth of smoking-induced mutations in 
> NSRO-driven non-small-cell lung cancer despite smoking exposure". 
> DOI: (for preprint) 10.1101/2023.07.04.547310

"NSRO" denotes non-smoking-related oncogene, the most prominent of which
is EGFR. In particular, KRAS is not in this group, because oncogenic
KRAS mutations are common in lung cancers in smokers.

The scripts are written in R (version 4.1.0). The data needed for the code
are in the file "ITH2_data.RData". 

The whole-exome and RNA sequencing data have been deposited at 
the European Genome-phenome Archive (EGA, http://www.ebi.ac.uk/ega/), 
under accession number EGAS00001006942.

# Links to external resources used in this study
* TCGA cBioPortal: http://www.cbioportal.org/
* Reactome pathways on MSigDB: https://www.gsea-msigdb.org/gsea/msigdb/human/genesets.jsp?collection=CP:REACTOME
* COSMIC catalog of driver gene mutations: https://cancer.sanger.ac.uk/cmc/home
* RNA gene panel used to classify gene expression subtype (TRU vs. non-TRU) by Wilkerson et al.: https://pubmed.ncbi.nlm.nih.gov/22590557/
