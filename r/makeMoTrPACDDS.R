# Edward Lau 2023
# Get MoTrPAC v8 all samples, save as DDS, then batch correct with ComBAT. This is for use with PLIER to harmonize MoTrPAC with GTEx

library(tidyverse)
library(sva)
library(DESeq2)
library(tibble)
library(org.Rn.eg.db)
library(SummarizedExperiment)

# Read the read count files for all four MoTrPAC tissues, then join into one data table

data_dir <- '../Data/MoTrPAC_Real/transcriptomics_results'
motr_samples <- c("T55_gastrocnemius", "T58_heart", "T68_liver", "T70_white_adipose")

# Open the metadata that will link vial ID to participant ID
metadata <- readr::read_tsv("../Data/MoTrPAC_Real/phenotype/merged_dmaqc_data2019-10-15.txt")
metadata$vial_label <- as.character(metadata$vial_label)
metadata$pid <- as.character(metadata$pid)
pidlist <- metadata %>% dplyr::pull(pid, name=vial_label)

# total dataframe
motr_all <- tibble()
subjects_all <- c()

for (sample_name in motr_samples) {
  short_sample_name <- gsub('T[0-9][0-9]_', '', sample_name)
  print(short_sample_name)
  motrpac <- readr::read_tsv(file.path(data_dir, sample_name, 'rna-seq', 'results', paste0('MoTrPAC_rsem_genes_count_', short_sample_name, '_v1.txt')))
  
  # Use AnnotationDbi to turn Ensembl Gene Name to Gene Symbol
  ens_to_gn = select(org.Rn.eg.db, keys=motrpac$gene_id, keytype='ENSEMBL', columns='SYMBOL') %>% dplyr::rename(gene_id = ENSEMBL)
  ens_to_gn <- ens_to_gn %>% dplyr::filter(!is.na(SYMBOL)) %>% dplyr::arrange(SYMBOL) %>% dplyr::distinct(gene_id, .keep_all=T)
  # Get gene name
  motrpac <- motrpac %>% dplyr::left_join(ens_to_gn)  %>% dplyr::filter(!is.na(SYMBOL)) %>% dplyr::select(-gene_id) %>% dplyr::relocate(SYMBOL)
  motrpac <- motrpac %>% dplyr::distinct(SYMBOL, .keep_all=T)
  # Transpose so we have n x p matrix
  
  # motrpac_t <- t(as.matrix(motrpac[,-1]))
  # colnames(motrpac_t) <- unlist(motrpac[, 1])
  # colnames(motrpac_t) <- paste0(short_sample_name, '_', colnames(motrpac_t))
  # 
  # # rownames(motrpac_t) <- pidlist[rownames(motrpac_t)]
  # 
  # motrpac_tibble <- as.tibble(motrpac_t)
  # 
  # # Get PID from Vial Label
  # motrpac_tibble$vial_label <- rownames(motrpac_t)
  # motrpac_tibble <- motrpac_tibble %>% dplyr::left_join(dplyr::select(metadata, vial_label, pid))
  # motrpac_tibble <- motrpac_tibble %>% dplyr::select(-vial_label) %>% dplyr::filter(!is.na(pid))
  #  
  # join to a full matrix by c_binding. First make sure all the matrices have the same rows
  if(nrow(motr_all) == 0){
    motr_all <- motrpac
  
  } else {
    motr_all <-  dplyr::left_join(motr_all, motrpac)
  }
  
}

# Get only samples that are in the study (no standards whose vial IDs start with 8)
motr_all <- dplyr::bind_cols(motr_all[1], motr_all[, colnames(motr_all) %in% metadata$vial_label])

# Make all gene names Upper case 
motr_all$SYMBOL <- toupper(motr_all$SYMBOL)

# Read the metadata, mainly we want to match the pid to animal agegroup and exercise group -- it appears the release 1 only has young age group for now (6 months)
table(metadata$animal.key.agegroup)
# For intervention, there are acute exercise (1) and control (3)
table(metadata$animal.key.intervention)
# For sacrifice time 
# |1|2|3|4|5|6|7|8|9|10|11
# |Immediate Post Exercise (Phase 1A)|0.5 hour (Phase 1A)|1 hour (Phase 1A)|4 hour (Phase 1A)|7 hour (Phase 1A)|24 hour (Phase 1A)|48 hour (Phase 1A)|1 week of training or control time (Phase 1B)|2 weeks of training (Phase 1B)|4 weeks of training (Phase 1B)|8 weeks of training or control time (Phase 1B)
table(metadata$animal.key.sacrificetime)

# There are two batch rows in the metadata dictionary, we can probably use registration.batchnumber for batch correction
table(metadata$animal.key.batch, metadata$animal.registration.batchnumber)

meta_df <- metadata %>% dplyr::select(vial_label, 
                                      pid, 
                                      animal.key.intervention, 
                                      animal.key.sacrificetime,
                                      animal.key.batch,
                                      animal.registration.batchnumber,
                                      animal.registration.sex,
                                      specimen.processing.sampletypedescription
                                      )

meta_df <- meta_df %>% column_to_rownames('vial_label')
meta_df$vial_label <- rownames(meta_df) 
meta_df <- meta_df %>% dplyr::relocate(vial_label)
meta_df <- meta_df[colnames(motr_all)[2:length(colnames(motr_all))], ]

# Turn into numerical matrix
motr_mat <- motr_all %>% tibble::column_to_rownames(var="SYMBOL") %>% as.matrix()
mode(motr_mat) <- "integer"

# Set up dds object
dds <- DESeq2::DESeqDataSetFromMatrix(countData = motr_mat,
                              colData = meta_df,
                              design = ~animal.registration.sex
                              )


saveRDS(dds, "../Projects/2022-04_PLIER/data/motrpac_dds.Rds")

# Normalize with vst
vsd <- DESeq2::vst(dds)
mat_vsd <- assay(vsd)

saveRDS(vsd, "../Projects/2022-04_PLIER/data/motrpac_vst_dds.Rds")

# https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2855-9
# I think we should first correct for technical batches, then correct for biological characteristics 
# 1. Animal experiment batch (animal.registration.batchnumber)
# 2. Sex (animal.registration.sex)


##
## Correct 1 - Batch correct with animal batch
##

mat_vsd_tec <- sva::ComBat(dat=mat_vsd, batch = meta_df$animal.registration.batchnumber) 

vsd_tec <- vsd
assay(vsd_tec) <- mat_vsd_tec
saveRDS(vsd_tec, "../Projects/2022-04_PLIER/data/motrpac_vst_combat1_dds.Rds")

##
## Correct 1 - Batch correct with extraction batch then sequencing batch
##
mat_vsd_tec_bio <- sva::ComBat(dat=mat_vsd_tec, batch = meta_df$animal.registration.sex) 

vsd_tec_bio <- vsd
assay(vsd_tec_bio) <- mat_vsd_tec_bio
saveRDS(vsd_tec_bio, "../Projects/2022-04_PLIER/data/motrpac_vst_combat2_dds.Rds")


