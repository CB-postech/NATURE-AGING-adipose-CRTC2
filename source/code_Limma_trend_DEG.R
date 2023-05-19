library(edgeR)
library(Seurat)
library(SingleCellExperiment)
library(scran)
library(scater)
setwd("/home/espark/project/mouse_CRTC/data/paper_final/")
seurat <- readRDS("seurat/qc_positive_all_seurat.rds")
#Subset seurat object by sample condition
seurat$figure_simple_annot <- as.factor(seurat$figure_simple_annot)
seurat_type <- SplitObject(seurat, split.by = "figure_simple_annot")
for (i in 1:nlevels(seurat$figure_simple_annot)){
  print(names(seurat_type[i]))
  sce <- as.SingleCellExperiment(seurat_type[[i]])
  #Normalizateion
  clusters <- quickCluster(sce)
  sce <- computeSumFactors(sce, clusters=clusters)
  summary(sizeFactors(sce))
  sce <- logNormCounts(sce)
  #Limma
  normcounts <- logcounts(sce)
  normcounts <- normcounts[rowSums(normcounts)!=0,]
  KO <- as.factor(colData(sce)[["KO"]])
  age <- as.factor(colData(sce)[["age"]])
  condition <-as.factor(colData(sce)[["condition"]])
  sample <- as.factor(colData(sce)[["sample"]])
  design <- model.matrix(~ 0+ condition)
  colnames(design) <- levels(condition)
  dupcor_remove <- duplicateCorrelation(normcounts, design, block=sample)
  fit_remove <- lmFit(normcounts, design, block = sample, correlation=dupcor_remove$consensus)
  cont.matrix <- makeContrasts(
    AGEdiff_inWT = Old_WT - Young_WT,
    AGEdiff_inKO = Old_KO - Young_KO,
    KOdiff_inYoung = Young_KO - Young_WT,
    KOdifF_inOld = Old_KO - Old_WT,
    Interaction = (Old_KO - Old_WT) - (Young_KO - Young_WT),
    levels = design
  )
  fit_remove2 <- contrasts.fit(fit_remove, cont.matrix)
  fit_remove2 <- eBayes(fit_remove2, trend = T, robust = T)
  saveRDS(fit_remove,file=paste0("limma_simple/",names(seurat_type[i]),'_limma_trend_result_scran.rds'))
  saveRDS(fit_remove2,file=paste0("limma/",names(seurat_type[i]),'_limma_trend_contrast_result_scran.rds'))
}

