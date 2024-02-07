# Run a Multiplier script
# This makes a PLIER model using the multiPLIER wrapper (this will be using a different set of priors and the recommended k value)

source("../Data/gtex_v8/plier/plier_util.R")

library(DESeq2)
library(dplyr)
library(readr)
library(tibble)
library(hpar)

data(hpaSecretome)

# Read the read count files
gtex <- readr::read_tsv("../Data/gtex_v8/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct", skip=2)

# Make distinct gene name
gtex_matrix <- gtex %>% dplyr::select(-Name) %>% dplyr::distinct(Description, .keep_all=TRUE) %>% tibble::column_to_rownames(var="Description") %>% as.matrix()

# Read the annotation file
annot <- readr::read_tsv("../Data/gtex_v8/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt", guess_max=15000)
annot$SUBJID <- gsub("(?<!GTEX)-.*$", "", annot$SAMPID, perl=T)

# Read the subject demographics file
subdemo <- readr::read_tsv("../Data/gtex_v8/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt")

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

# vsd <- DESeq2::vst(dds)

gtex_data <- assay(dds)

# Remove Secretome
hpa_secreted_to_blood <- hpaSecretome %>% dplyr::filter(`Secretome.location` == "Secreted to blood") %>% dplyr::select(`Gene.name`, `Secretome.location`)
gtex_data_sans_secretome <- gtex_data[!rownames(gtex_data) %in% hpa_secreted_to_blood$Gene.name, ]



plierResults_gtexVst_bloodCell_svm_canonicalPathways <- PLIERNewData(gtex_data_sans_secretome)

saveRDS(plierResults_gtexVst_bloodCell_svm_canonicalPathways, "../Data/gtex_v8/plier/GTExv8RawNoSecretome_canonicalPathwaysBloodCellSVM_multiplier.Rds")