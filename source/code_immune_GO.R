library(edgeR)
library(Seurat)
myeloid_seurat <- readRDS("data/seurat/qc_positive_myeloid_seurat.rds")
lymphoid_seurat <- readRDS("data/seurat/qc_positive_lymphoid_seurat.rds")
lymphoid_seurat <- subset(lymphoid_seurat, cells=colnames(lymphoid_seurat)[
  lymphoid_seurat$lymphoid_annot!="T.prolif"])


#Load Lymphoid cell type limma trend result
df_all_age <- data.frame()
df_gene_age <- data.frame()
for (i in 1:nlevels(myeloid_seurat$myeloid_annot)){
  print(levels(myeloid_seurat$myeloid_annot)[i])
  assign("tmp",readRDS(paste0("limma/",levels(myeloid_seurat$myeloid_annot)[i],"_limma_trend_contrast_result_scran.rds")))
  assign(paste0(levels (myeloid_seurat$myeloid_annot)[i],"_result"),tmp)
  #Difference Old Young WT DEG
  assign("tmp1", topTable(tmp, coef = 1, sort.by = "P", n = Inf))
  assign(paste0("table_",levels(myeloid_seurat$myeloid_annot)[i]),tmp1)
  tmp1_df <- data.frame(genes=rownames(tmp1), logFC=tmp1$logFC, cell_type= levels(myeloid_seurat$myeloid_annot)[i], p_value =tmp1$adj.P.Val)
  df_all_age <- rbind(df_all_age,tmp1_df)
  head(df_all_age)
  
  #Extract only significant DEG
  assign("tmp2", data.frame(tmp1)[tmp1$adj.P.Val<0.05,])
  if( length(rownames(tmp2)) != 0)
  {
    tmp_df <-data.frame(genes=rownames(tmp2),cell_type=levels(myeloid_seurat$myeloid_annot)[i], log_FC= tmp2$logFC,
                        adj.P.Val = tmp2$adj.P.Val)
    df_gene_age <- rbind(df_gene_age,tmp_df)  
  }
}
for (i in 1:nlevels(lymphoid_seurat$lymphoid_annot)){
  print(levels(lymphoid_seurat$lymphoid_annot)[i])
  assign("tmp",readRDS(paste0("limma/",levels(lymphoid_seurat$lymphoid_annot)[i],"_limma_trend_contrast_result_scran.rds")))
  assign(paste0(levels (lymphoid_seurat$lymphoid_annot)[i],"_result"),tmp)
  #Difference Old Young WT DEG
  assign("tmp1", topTable(tmp, coef = 1, sort.by = "P", n = Inf))
  assign(paste0("table_",levels(lymphoid_seurat$lymphoid_annot)[i]),tmp1)
  tmp1_df <- data.frame(genes=rownames(tmp1), logFC=tmp1$logFC, cell_type= levels(lymphoid_seurat$lymphoid_annot)[i], p_value =tmp1$adj.P.Val)
  df_all_age <- rbind(df_all_age,tmp1_df)
  head(df_all_age)
  
  #Extract only significant DEG
  assign("tmp2", data.frame(tmp1)[tmp1$adj.P.Val<0.05,])
  if( length(rownames(tmp2)) != 0)
  {
    tmp_df <-data.frame(genes=rownames(tmp2),cell_type=levels(lymphoid_seurat$lymphoid_annot)[i], log_FC= tmp2$logFC,
                        adj.P.Val = tmp2$adj.P.Val)
    df_gene_age <- rbind(df_gene_age,tmp_df)  
  }
}
up_DEG_WT <- df_gene_age[df_gene_age$log_FC>0,]
down_DEG_WT <-df_gene_age[df_gene_age$log_FC<0,]

up_GO_WT <- data.frame()
down_GO_WT <- data.frame()
i= levels(up_DEG_WT$cell_type)[1]
for (i in levels(up_DEG_WT$cell_type)){
  print(i)
  tmp_DEG <- up_DEG_WT[up_DEG_WT$cell_type==i,]
  tmp_GO <- topGO_function(tmp_DEG$genes, rownames(myeloid_seurat))
  tmp_GO$celltype <- i
  up_GO_WT <- rbind(up_GO_WT,tmp_GO)
}
for (i in levels(down_DEG_WT$cell_type)[c(10:14,16:19)]){
  print(i)
  tmp_DEG <- down_DEG_WT[down_DEG_WT$cell_type==i,]
  tmp_GO <- topGO_function(tmp_DEG$genes, rownames(myeloid_seurat))
  tmp_GO$celltype <- i
  down_GO_WT <- rbind(down_GO_WT,tmp_GO)
}
