# Run a Multiplier script
# This makes a PLIER model using the multiPLIER wrapper (this will be using a different set of priors and the recommended k value)

source("../Data/gtex_v8/plier/plier_util.R")
library(DESeq2)

gtex_vst_combat <- readRDS("../Data/gtex_v8/gtexv8_vst_combat3_dds.Rds")

plierResults_gtexVstCombat_bloodCell_svm_canonicalPathways <- PLIERNewData(assay(gtex_vst_combat))

saveRDS(plierResults_gtexVstCombat_bloodCell_svm_canonicalPathways, "../Data/gtex_v8/plier/GTExv8VSTcomBAT_canonicalPathwaysBloodCellSVM_multiplier.Rds")