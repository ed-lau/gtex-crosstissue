# Run a Multiplier script
# 2023-02-28 Age corrected
# This makes a PLIER model using the multiPLIER wrapper (this will be using a different set of priors and the recommended k value)

source("../Data/gtex_v8/plier/plier_util.R")

library(DESeq2)
library(dplyr)
library(readr)
library(tibble)
library(hpar)

data(hpaSecretome)

gtex_vst_combat <- readRDS("../Data/gtex_v8/gtexv8_vst_combat4_dds.Rds")
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

# Overwrite PLIER function because the way it checkes class doesn't work with R 4.0.0 matrix objects
#' estimates the number of 'significant' principle components for the SVD decomposition -- this is the minimum k for PLIER

#' @param  data the same data as to be used for PLIER (z-score recommended) or alternatively the result of an svd calculation 
#' @param method Either "eblow" (fast) or "permutation" (slower, but less heuristic)
#' @param B number of permutations
#' @param seed seed for reproducibility 
#' @export
num.pc2 = function (data, method="elbow", B = 20, seed = NULL) 
{
  
  method=match.arg(method, c("elbow", "permutation"))
  if (!is.null(seed)) {
    set.seed(seed)
  }
  warn <- NULL
  if(TRUE){
  #if((class(data)!="list")&(class(data)!="rsvd")){
    message("Computing svd")
    n <- ncol(data)
    m <- nrow(data)
    data=rowNorm(data)
    if(n<500){
      k=n
    }
    else{
      k=max(200,n/4)
    }
    if(k==n){
      uu <- svd(data)
    }
    else{
      set.seed(123456);uu <- rsvd(data,k, q=3)
    }
  }
  else if (!is.null(data[["d"]])){
    if(method=="permutation"){
      message("Original data is needed for permutation method.\nSetting method to elbow")
      method="elbow"
    }
    
    uu=data
    
  }
  
  
  
  if(method=="permutation"){
    #nn = min(c(n, m))
    dstat <- uu$d[1:k]^2/sum(uu$d[1:k]^2)
    dstat0 <- matrix(0, nrow = B, ncol = k)
    for (i in 1:B) {
      dat0 <- t(apply(data, 1, sample, replace = FALSE))
      if(k==n){
        uu0 <- svd(dat0)
      }
      else{
        set.seed(123456);
        uu0 <- rsvd(dat0,k, q=3)
      }
      dstat0[i, ] <- uu0$d[1:k]^2/sum(uu0$d[1:k]^2)
    }
    psv <- rep(1, k)
    for (i in 1:k) {
      psv[i] <- mean(dstat0[, i] >= dstat[i])
    }
    for (i in 2:k) {
      psv[i] <- max(psv[(i - 1)], psv[i])
    }
    
    nsv <- sum(psv <= 0.1)
  }
  else if (method=="elbow"){
    x=smooth(xraw<-abs(diff(diff(uu$d))), twiceit = T)
    #plot(x)
    
    
    nsv=which(x<=quantile(x, 0.5))[1]+1
    
  }
  return(nsv)
}

PLIERNewData2 <- function(exprs.mat, seed = 12345) {
  # A wrapper function for applying PLIER to a data set. We use the following
  # genesets that come with PLIER: bloodCellMarkersIRISDMAP, svmMarkers, 
  # and canonicalPathways. We set the k parameter for the PLIER model by
  # identifying the number of "significant PCs" with PLIER::num.pc and then 
  # using sig PCs * 0.3. This is consistent with recommendations from the 
  # PLIER authors.
  # 
  # Args:
  #   exprs.mat: a gene expression matrix, rows are genes, columns are samples
  #   seed: an integer to be supplied to set.seed() for reproducibility 
  #         purposes, default is 12345
  #         
  # Returns:
  #   plier.res: output from PLIER::PLIER()
  #
  require(PLIER)
  
  set.seed(seed)
  
  # load PLIER pathway and cell type data
  data(bloodCellMarkersIRISDMAP)
  data(svmMarkers)
  data(canonicalPathways)
  
  # combine the pathway data from PLIER
  all.paths <- PLIER::combinePaths(bloodCellMarkersIRISDMAP, svmMarkers, 
                                   canonicalPathways)
  
  # what genes are common to the pathway data and the expression matrix
  cm.genes <- PLIER::commonRows(all.paths, exprs.mat)
  
  # row normalize
  exprs.norm <- PLIER::rowNorm(exprs.mat)
  exprs.data <- exprs.norm[cm.genes, ]

  # what should we set the minimum k parameter to in PLIER? estimate the number 
  # of PC for the SVD decomposition 
  set.k <- num.pc2(exprs.data)
  
  # PLIER main function + return results
  plier.res <- PLIER::PLIER(exprs.norm[cm.genes, ], all.paths[cm.genes, ], 
                            k = round((set.k + set.k * 0.3), 0), trace = TRUE)
  
  return(plier.res)
  
}

# Forking our own version of the num.pc function in PLIER because it is giving a condition for length > 1 error
plierResults_gtexVstCombatNoSec_bloodCell_svm_canonicalPathways <- PLIERNewData2(gtex_vst_comBAT_nosecretome)

saveRDS(plierResults_gtexVstCombatNoSec_bloodCell_svm_canonicalPathways, "../Data/gtex_v8/plier/GTExv8VSTcomBATAgeCorrNoSecretome_canonicalPathwaysBloodCellSVM_multiplier_v2
        .Rds")