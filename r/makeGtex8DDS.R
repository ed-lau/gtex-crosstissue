# Edward Lau 2021
# Get GTEx v8 all samples, save as DDS, then batch correct with ComBAT.

library(tidyverse)
library(sva)
library(DESeq2)


# Read the read count files
gtex <- readr::read_tsv("../Edward/gtex_v8/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct", skip=2)

# Make distinct gene name
gtex_matrix <- gtex %>% dplyr::select(-Name) %>% dplyr::distinct(Description, .keep_all=TRUE) %>% tibble::column_to_rownames(var="Description") %>% as.matrix()

# Read the annotation file
annot <- readr::read_tsv("../Edward/gtex_v8/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt", guess_max=15000)
annot$SUBJID <- gsub("(?<!GTEX)-.*$", "", annot$SAMPID, perl=T)

# Read the subject demographics file
subdemo <- readr::read_tsv("../Edward/gtex_v8/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt")

# Get the metadata file for all RNA-seq data
i. <- annot %>% dplyr::filter(SMAFRZE == "RNASEQ")
i. <- i. %>% dplyr::left_join(subdemo)
i.$rowname <- i.$SAMPID
i. <- i. %>% column_to_rownames("rowname")
i. <- i.[colnames(gtex_matrix),]

# Set up dds object
dds <- DESeqDataSetFromMatrix(countData = gtex_matrix,
                              colData = i.,
                              design= ~ SMTSD)

# Normalize with vst
vsd <- DESeq2::vst(dds)
mat_vsd <- assay(vsd)

saveRDS(vsd, "../Edward/gtex_v8/gtexv8_vst_dds.Rds")

# https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2855-9
# I think we should first correct for technical batches, then correct for biological characteristics 
# 1. Extraction Batch (SMNABTCH)
# 2. Sequencing Batch (SMGEBTCH)
# 3. Tissue Ischemic time (bin into 300 min intervals) (SMTSISCH) cut(i..$SMTSISCH, breaks=c(-Inf, 0, 300, 600, 900, 1200, Inf)) %>% table()
# 4. Death type (violent and quick vs. slow) (DTHHRDY) 
# 5. 2023-02-28, add AGE for crosstalk.

##
## Correct 1 - Batch correct with extraction batch then sequencing batch
##
correct_technical <- i. %>% dplyr::select(SUBJID, SMNABTCH, SMGEBTCH)
correct_technical <- correct_technical[colnames(mat_vsd),]

        # table(correct_technical$SMNABTCH)
        # table(correct_technical$SMGEBTCH)

correct_technical$SMGEBTCH[is.na(correct_technical$SMGEBTCH)] <- "NoBatch"
mat_vsd_tec <- sva::ComBat(dat=mat_vsd, batch = correct_technical$SMNABTCH) 
mat_vsd_tec <- sva::ComBat(dat=mat_vsd_tec, batch = correct_technical$SMGEBTCH) 

vsd_tec <- vsd
assay(vsd_tec) <- mat_vsd_tec
saveRDS(vsd_tec, "../Edward/gtex_v8/gtexv8_vst_combat1_dds.Rds")

# Histogram of VSD
#h._vsd_tech[,1] %>% hist()
# Plot of pre-and-post correction
#plot(h._vsd_tech[,1], h._vsd[,1])

##
## Correct 2  - Batch correct with Tissue biological characteristics (ischemia time, potentially add fixative and autolysis score)
##
correct_tissue <- i. %>% dplyr::select(SUBJID, SMTSISCH, SMTSPAX, SMATSSCR, SMTSD)
correct_tissue <- correct_tissue[colnames(mat_vsd),]
correct_tissue$ischemic_batch <- cut(correct_tissue$SMTSISCH, breaks=c(-Inf, 300, 600, 900, 1200, Inf)) 
#correct_tissue$fixative_batch <- cut(correct_tissue$SMTSPAX, breaks=c(300, 600, 900, 1200, Inf)) 
#correct_tissue$autolysis_batch <- cut(correct_tissue$SMATSSCR, breaks=c(300, 600, 900, 1200, Inf)) 

# Check if ischemic time correlated with organ
        # ggplot(data=correct_tissue, aes(x=SMTSD, y=SMTSISCH)) + geom_boxplot()

# Check that we have a sensible distribution of the batches
        # table(correct_tissue$ischemic_batch)

# If no information, encode as no specific ischemia time.
correct_tissue$ischemic_batch[is.na(correct_tissue$ischemic_batch)] <- "(-Inf,300]"

mat_vsd_tec_tis <- sva::ComBat(dat=mat_vsd_tec, batch = correct_tissue$ischemic_batch) 

# Histogram of VSD
        # mat_vsd_tec_tis[,1] %>% hist()
# Plot of pre-and-post correction
        # plot(h._vsd_tech[,1], h._vsd_tech_tiss[,1])

vsd_tec_tis <- vsd
assay(vsd_tec_tis) <- mat_vsd_tec_tis
saveRDS(vsd_tec_tis, "../Edward/gtex_v8/gtexv8_vst_combat2_dds.Rds")


# Correct 3 - Batch correct with donor death type
correct_donor <- tibble::enframe(colnames(mat_vsd)) %>% dplyr::select(SAMPID = value)
correct_donor$SUBJID <- gsub("(?<!GTEX)-.*$", "", correct_donor$SAMPID, perl=T)
correct_donor <- correct_donor %>% dplyr::left_join(subdemo) %>% tibble::column_to_rownames("SAMPID")
correct_donor <- correct_donor[colnames(mat_vsd),]
# Check that we have a sensible distribution of the batches
        # table(correct_donor$DTHHRDY)
        # table(correct_donor$AGE)
correct_donor$DTHHRDY[is.na(correct_donor$DTHHRDY)] <- 0
mat_vsd_tec_tis_bio <- sva::ComBat(dat=mat_vsd_tec_tis, batch = correct_donor$DTHHRDY)

# Histogram of VSD
        # mat_vsd_tec_tis_bio[,1] %>% hist()
# Plot of pre-and-post correction
        # plot( mat_vsd_tec_tis[,1], mat_vsd_tec_tis_bio[,1])

# Check whether there are genes with almost 0 variances (disable comment to run since it takes a long time, but last I check there was no nearZeroVar genes)
        # dim(mat_vsd_tec_tis_bio)
        # caret::nearZeroVar(t(h._vsd_tech_tiss_biol))

vsd_tec_tis_bio <- vsd
assay(vsd_tec_tis_bio) <- mat_vsd_tec_tis_bio
saveRDS(vsd_tec_tis_bio, "../Edward/gtex_v8/gtexv8_vst_combat3_dds.Rds")


# NEW 2023-02-28 - We are going to try to correct for age for the R03 tissue crosstalk project
correct_donor$AgeGroup <- as.integer(gsub("-.*$", "", correct_donor$AGE, ))
mat_vsd_tec_tis_bio_age <- sva::ComBat(dat=mat_vsd_tec_tis_bio, batch = correct_donor$AgeGroup)

vsd_tec_tis_bio_age <- vsd
assay(vsd_tec_tis_bio_age) <- mat_vsd_tec_tis_bio_age
saveRDS(vsd_tec_tis_bio_age, "../Edward/gtex_v8/gtexv8_vst_combat4_dds.Rds")



