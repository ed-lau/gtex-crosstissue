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

# Read the read count files and make the uncorrected GTEx matrix
h <- readr::read_tsv("/../gtex_v8/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct", skip=2)
gtex_matrix <- h %>% dplyr::select(-Name) %>% dplyr::distinct(Description, .keep_all=TRUE) %>% tibble::column_to_rownames(var="Description") %>% as.matrix()

# Read the annotation file
annot <- readr::read_tsv("../gtex_v8/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt", guess_max=15000)
annot$SUBJID <- gsub("(?<!GTEX)-.*$", "", annot$SAMPID, perl=T)

# Read the subject demographics file
subdemo <- readr::read_tsv("../Data/gtex_v8/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt")

# Get the metadata file for all RNA-seq data
i. <- annot %>% dplyr::filter(SMAFRZE == "RNASEQ")
i. <- i. %>% dplyr::left_join(subdemo)
i.$rowname <- i.$SAMPID
i. <- i. %>% column_to_rownames("rowname")
i. <- i.[colnames(gtex_matrix),]

gtex <- DESeqDataSetFromMatrix(countData = gtex_matrix,
                              colData = i.,
                              design= ~ SMTSD)

# Read the saved VST and ComBAT objects
vsd <- readRDS("../gtex_v8/gtexv8_vst_dds.Rds")
vsd_tec <- readRDS("../gtex_v8/gtexv8_vst_combat1_dds.Rds")
vsd_tec_tis <- readRDS("../gtex_v8/gtexv8_vst_combat2_dds.Rds")
vsd_tec_bio <- readRDS("../gtex_v8/gtexv8_vst_combat3_dds.Rds")
vsd_tec_bio_age <- readRDS("../gtex_v8/gtexv8_vst_combat4_dds.Rds")

# We should do a PCA plot before and after each step of normalization/ batch correction
p0 <- rsvd::rpca(assay(gtex), k=2)
p0_df <- p0$rotation %>% as.data.frame() %>% tibble::rownames_to_column() %>% tibble::as_tibble() %>% dplyr::select(SAMPID = rowname, PC1:PC2) %>% dplyr::left_join(i.)
g0 <- ggplot(data=p0_df, aes(x=PC1, y=PC2)) 
g0 <- g0 + geom_point(aes(col=SMTS), size=0.7, stroke=0) + ggthemes::theme_few() + theme(aspect.ratio=1, legend.position="bottom")
g0 <- g0 + ggtitle("Raw Counts")

p1 <- rsvd::rpca(assay(vsd), k=2)
p1_df <- p1$rotation %>% as.data.frame() %>% tibble::rownames_to_column() %>% tibble::as_tibble() %>% dplyr::select(SAMPID = rowname, PC1:PC2) %>% dplyr::left_join(i.)
g1 <- ggplot(data=p1_df, aes(x=PC1, y=PC2))
g1 <- g1 + geom_point(aes(col=SMTS), size=0.7, stroke=0) + ggthemes::theme_few() + theme(aspect.ratio=1, legend.position="bottom")
g1 <- g1 + ggtitle("After Variance Stabilization")

p2 <- rsvd::rpca(assay(vsd_tec), k=2)
p2_df <- p2$rotation %>% as.data.frame() %>% tibble::rownames_to_column()  %>% tibble::as_tibble() %>% dplyr::select(SAMPID = rowname, PC1:PC2) %>% dplyr::left_join(i.)
g2 <- ggplot(data=p2_df, aes(x=PC1, y=PC2))
g2 <- g2 + geom_point(aes(col=SMTS), size=0.7, stroke=0) + ggthemes::theme_few() + theme(aspect.ratio=1, legend.position="bottom")
g2 <- g2 + ggtitle("After Stepwise Technical Batch Correction")

p3 <- rsvd::rpca(assay(vsd_tec_tis), k=2)
p3_df <- p3$rotation %>% as.data.frame() %>% tibble::rownames_to_column()  %>% tibble::as_tibble() %>% dplyr::select(SAMPID = rowname, PC1:PC2) %>% dplyr::left_join(i.)
g3 <- ggplot(data=p3_df, aes(x=PC1, y=PC2))
g3 <- g3 + geom_point(aes(col=SMTS), size=0.7, stroke=0) + ggthemes::theme_few() + theme(aspect.ratio=1, legend.position="bottom")
g3 <- g3 + ggtitle("After Stepwise Tissue Charateristics Correction")

p4 <- rsvd::rpca(assay(vsd_tec_bio), k=2)
p4_df <- p4$rotation %>% as.data.frame() %>% tibble::rownames_to_column()  %>% tibble::as_tibble() %>% dplyr::select(SAMPID = rowname, PC1:PC2) %>% dplyr::left_join(i.)
g4 <- ggplot(data=p4_df, aes(x=PC1, y=PC2))
g4 <- g4 + geom_point(aes(col=SMTS), size=0.7, stroke=0) + ggthemes::theme_few() + theme(aspect.ratio=1, legend.position="bottom")
g4 <- g4 + ggtitle("After Stepwise Donor Death Charateristics Correction")

p5 <- rsvd::rpca(assay(vsd_tec_bio_age), k=2)
p5_df <- p5$rotation %>% as.data.frame() %>% tibble::rownames_to_column()  %>% tibble::as_tibble() %>% dplyr::select(SAMPID = rowname, PC1:PC2) %>% dplyr::left_join(i.)
g5 <- ggplot(data=p5_df, aes(x=PC1, y=PC2))
g5 <- g5 + geom_point(aes(col=SMTS), size=0.7, stroke=0) + ggthemes::theme_few() + theme(aspect.ratio=1, legend.position="bottom")
g5 <- g5 + ggtitle("After Stepwise Age Correction")

ggsave(plot=g0, file= "../figs/GTExv8_batchCorrection_1.pdf", width=5, height=5, useDingbats=F)
ggsave(plot=g1, file= "../figs/GTExv8_batchCorrection_2.pdf", width=5, height=5, useDingbats=F)
ggsave(plot=g2, file= "../figs/GTExv8_batchCorrection_3.pdf", width=5, height=5, useDingbats=F)
ggsave(plot=g3, file= "../figs/GTExv8_batchCorrection_4.pdf", width=5, height=5, useDingbats=F)
ggsave(plot=g4, file= "../figs/GTExv8_batchCorrection_5.pdf", width=5, height=5, useDingbats=F)
ggsave(plot=g5, file= "../figs/GTExv8_batchCorrection_6.pdf", width=5, height=5, useDingbats=F)
