# Run a Multiplier script
# This makes a PLIER model using the multiPLIER wrapper (this will be using a different set of priors and the recommended k value)

source("../Data/gtex_v8/plier/plier_util.R")

library(DESeq2)
library(dplyr)
library(readr)
library(tibble)
library(hpar)

data(hpaSecretome)

gtex_vst_combat <- readRDS("../Data/gtex_v8/gtexv8_vst_combat3_dds.Rds")
gtex_vst_combat <- assay(gtex_vst_combat)


# ligand_annotations = OmnipathR::import_omnipath_annotations(resources=c('HPA_secretome', 'UniProt_keyword'))
# ligand_annotations_HPAligands <- ligand_annotations %>% dplyr::filter(source == 'HPA_secretome', entity_type == "protein", value == "Secreted to blood")
# ligand_annotations_UniProt <- ligand_annotations %>% dplyr::filter(source == 'UniProt_keyword', entity_type == "protein", value=="Secreted")

# Combine HPA secretome and UniProt secreted to build a common list (Note: but maybe this list is too long and has about 2000 genes, maybe
# we should use a mroe stringent list for now)
# all_secreted <- unique(c(ligand_annotations_HPAligands$genesymbol, ligand_annotations_UniProt$genesymbol))


# commonGenes <- PLIER::commonRows(ligand_annotations_HPAligands %>% dplyr::distinct(genesymbol, .keep_all=T) %>% tibble::column_to_rownames("genesymbol"), gtex_vst_combat)
# non_secretome <- rownames(gtex_vst_combat)[!rownames(gtex_vst_combat) %in% commonGenes]
# Total 695 high-stringency endocrine candidates


# Remove Secretome
hpa_secreted_to_blood <- hpaSecretome %>% dplyr::filter(`Secretome.location` == "Secreted to blood") %>% dplyr::select(`Gene.name`, `Secretome.location`)
gtex_vst_comBAT_nosecretome <- gtex_vst_combat[!rownames(gtex_vst_combat) %in% hpa_secreted_to_blood$Gene.name, ]



plierResults_gtexVstCombatNoSec_bloodCell_svm_canonicalPathways <- PLIERNewData(gtex_vst_comBAT_nosecretome)

saveRDS(plierResults_gtexVstCombatNoSec_bloodCell_svm_canonicalPathways, "../Data/gtex_v8/plier/GTExv8VSTcomBATNoSecretome_canonicalPathwaysBloodCellSVM_multiplier_v2
        .Rds")