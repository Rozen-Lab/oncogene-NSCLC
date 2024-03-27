# About this code repository
This repository provides R codes and data used for the paper

> Huang, C.-Y., et al. (2023). "Dearth of Smoking-induced Mutations in Oncogene-Driven NSCLCs Despite Smoking Exposure". 
> DOI: (for preprint) 10.1101/2023.07.04.547310

The scripts are written in R (version 4.1.0). The data needed for the code
are in the file "ITH2_data.RData". 

The whole-exome and RNA sequencing data have been deposited at 
the European Genome-phenome Archive (EGA, http://www.ebi.ac.uk/ega/), 
under accession number EGAS00001006942.

# Links to external resources used in this study
* TCGA-LUAD: US National Cancer Institute GDC (Genomic Data Commons) Data Portal at https://portal.gdc.cancer.gov/projects/TCGA-LUAD (dbGap Study Accession phs000178, https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000178.v11.p8) cBioPortal: http://www.cbioportal.org/
* Chen et al., 2020, doi: 10.1038/s41588-019-0569-6, sequencing data: European Genome-Phenome Archive, accession numbers GAD00001004421 and EGAD00001004422 at https://ega-archive.org/studies/EGAS00001002941
* Chen et al., 2020, doi: 10.1038/s41588-019-0569-6, somatic mutations, clinical information: the OncoSG Cancer Genomics Portal, accession “Lung Adenocarcinoma (GIS, 2019)”, https://src.gisapps.org/OncoSG_public/study/summary?id=GIS031
* Reactome pathways on MSigDB: https://www.gsea-msigdb.org/gsea/msigdb/human/genesets.jsp?collection=CP:REACTOME
* COSMIC catalog of driver gene mutations: https://cancer.sanger.ac.uk/cmc/home
* RNA gene panel used to classify gene expression subtype (TRU vs. non-TRU) by Wilkerson et al.: https://pubmed.ncbi.nlm.nih.gov/22590557/
