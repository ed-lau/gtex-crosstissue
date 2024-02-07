# Edward Lau 2021
# Plot out the batch effects of GTEx v8 before and after ComBAT using fast rpca.

library(dplyr)
library(tibble)
library(readr)
library(ggplot2)
library(ggthemes)
library(DESeq2)
library(rsvd)
library(ggpubr)


# Read the saved DDS VST and ComBAT objects
dds <- readRDS("../Projects/2022-04_PLIER/data/motrpac_dds.Rds")
vsd <- readRDS("../Projects/2022-04_PLIER/data/motrpac_vst_dds.Rds")
vsd_tec <- readRDS("../Projects/2022-04_PLIER/data/motrpac_vst_combat1_dds.Rds")
vsd_tec_bio <- readRDS("../Projects/2022-04_PLIER/data/motrpac_vst_combat2_dds.Rds")

metadata  <- readr::read_tsv("../Data/MoTrPAC_Real/phenotype/merged_dmaqc_data2019-10-15.txt")
metadata$vial_label <- as.character(metadata$vial_label)

# We should do a PCA plot before and after each step of normalization/ batch correction
p0 <- rsvd::rpca(assay(dds), k=2)
p0_df <- p0$rotation %>% as.data.frame() %>% tibble::rownames_to_column() %>% tibble::as_tibble() %>% dplyr::select(vial_label = rowname, PC1:PC2) %>% dplyr::left_join(metadata)
g0 <- ggplot(data=p0_df, aes(x=PC1, y=PC2)) 
g0 <- g0 + geom_point(aes(col=specimen.processing.sampletypedescription), size=0.7, stroke=0) + ggthemes::theme_few() + theme(aspect.ratio=1, legend.position="bottom")
g0 <- g0 + ggtitle("Raw Counts")

p1 <- rsvd::rpca(assay(vsd), k=2)
p1_df <- p1$rotation %>% as.data.frame() %>% tibble::rownames_to_column() %>% tibble::as_tibble() %>% dplyr::select(vial_label = rowname, PC1:PC2) %>% dplyr::left_join(metadata)
g1 <- ggplot(data=p1_df, aes(x=PC1, y=PC2))
g1 <- g1 + geom_point(aes(col=specimen.processing.sampletypedescription), size=0.7, stroke=0) + ggthemes::theme_few() + theme(aspect.ratio=1, legend.position="bottom")
g1 <- g1 + ggtitle("After Variance Stabilization")

p2 <- rsvd::rpca(assay(vsd_tec), k=2)
p2_df <- p2$rotation %>% as.data.frame() %>% tibble::rownames_to_column() %>% tibble::as_tibble() %>% dplyr::select(vial_label = rowname, PC1:PC2) %>% dplyr::left_join(metadata)
g2 <- ggplot(data=p2_df, aes(x=PC1, y=PC2))
g2 <- g2 + geom_point(aes(col=specimen.processing.sampletypedescription), size=0.7, stroke=0) + ggthemes::theme_few() + theme(aspect.ratio=1, legend.position="bottom")
g2 <- g2 + ggtitle("After Stepwise Animal Batch Correction")

p4 <- rsvd::rpca(assay(vsd_tec_bio), k=2)
p4_df <- p4$rotation %>% as.data.frame() %>% tibble::rownames_to_column() %>% tibble::as_tibble() %>% dplyr::select(vial_label = rowname, PC1:PC2) %>% dplyr::left_join(metadata)
g4 <- ggplot(data=p4_df, aes(x=PC1, y=PC2))
g4 <- g4 + geom_point(aes(col=specimen.processing.sampletypedescription), size=0.7, stroke=0) + ggthemes::theme_few() + theme(aspect.ratio=1, legend.position="bottom")
g4 <- g4 + ggtitle("After Stepwise Sex Correction")

# p <- ggpubr::ggarrange(g0, g1, g2, g4, ncol=2, common.legend = T, legend="bottom")

ggsave(plot=g0, file= "../Projects/2022-04_PLIER/figs/MoTrPAC_batchCorrection_1.pdf", width=5, height=5, useDingbats=F)
ggsave(plot=g1, file= "../Projects/2022-04_PLIER/figs/MoTrPAC_batchCorrection_2.pdf", width=5, height=5, useDingbats=F)
ggsave(plot=g2, file= "../Projects/2022-04_PLIER/figs/MoTrPAC_batchCorrection_3.pdf", width=5, height=5, useDingbats=F)
ggsave(plot=g4, file= "../Projects/2022-04_PLIER/figs/MoTrPAC_batchCorrection_4.pdf", width=5, height=5, useDingbats=F)
