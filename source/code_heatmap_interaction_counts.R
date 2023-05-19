library(Seurat)
Young_WT_count_matrix <- read.csv("cellphoneDB/Young_WT_count_matrix.csv",row.names = 1)
Young_KO_count_matrix <- read.csv("cellphoneDB/Young_KO_count_matrix.csv",row.names = 1)
Old_WT_count_matrix <- read.csv("cellphoneDB/Old_WT_count_matrix.csv",row.names = 1)
Old_KO_count_matrix <-read.csv("cellphoneDB/Old_KO_count_matrix.csv",row.names = 1)
qc_crtc_seurat <- readRDS("data/seurat/qc_positive_all_seurat.rds")
# qc_crtc_seurat <- APC_seurat
max_count <- pmax(max(Young_WT_count_matrix),max(Young_KO_count_matrix),max(Old_WT_count_matrix),max(Old_KO_count_matrix))
library(RColorBrewer)
palette_length = 100
my_color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(palette_length)

my_breaks <- c(seq(0,max_count,
                   length.out=ceiling(palette_length)))

library(dplyr)
age_diff_count_matrix = Old_WT_count_matrix - Young_WT_count_matrix
KO_diff_count_matrix =  Old_KO_count_matrix- Old_WT_count_matrix

age_diff_count_matrix <- age_diff_count_matrix[-1,-1]
KO_diff_count_matrix <- KO_diff_count_matrix[-1,-1]


library(RColorBrewer)
palette_length = 100
my_color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(palette_length)
my_breaks <- c(seq(-27, 0,
                   length.out=ceiling(palette_length/2) + 1),
               seq(33/palette_length,
                   33,
                   length.out=floor(palette_length/2)))
g <-pheatmap::pheatmap(age_diff_count_matrix[], color = my_color,cluster_rows = F,cluster_cols = F,
                       breaks = my_breaks,border_color = NA)

gene_annot_df <- data.frame(row.names = c(APC,myeloid,lymphoid))
gene_annot_df$annot <- c(rep("Nonimmune",length(APC)),rep("Myeloid", length(myeloid)),
                         rep("Lymphoid", length(lymphoid)))
annot_col_df <- data.frame(row.names = c("Nonimmune","Myeloid","Lymphoid"))
annot_col_df <- c("#487947","#8A4512","#0070C0")
annot_col <- list(annot=c(Nonimmune="#8A4512",Myeloid="#487947",Lymphoid="#0070C0"))
check <- age_diff_count_matrix[c(APC,myeloid,lymphoid),c(APC,myeloid,lymphoid)]
for( i in 1:nrow(check)){
  print(i)
  for (j in 1:ncol(check)){
    if(i>j){
      check[i,j] <-NA
    }
  }
}
g <-pheatmap::pheatmap(age_diff_count_matrix[c(APC,myeloid,lymphoid),c(APC,myeloid,lymphoid)], color = my_color,cluster_rows = F,cluster_cols = F,
                       breaks = my_breaks,border_color = NA,annotation_row = gene_annot_df,
                       annotation_col = gene_annot_df,
                       annotation_colors = annot_col)


check <- KO_diff_count_matrix[c(APC,myeloid,lymphoid),c(APC,myeloid,lymphoid)]
for( i in 1:nrow(check)){
  print(i)
  for (j in 1:ncol(check)){
    if(i>j){
      check[i,j] <-NA
    }
  }
}
g <-pheatmap::pheatmap(KO_diff_count_matrix[c(APC,myeloid,lymphoid),c(APC,myeloid,lymphoid)], color = my_color,cluster_rows = F,cluster_cols = F,
                       breaks = my_breaks,border_color = NA,annotation_row = gene_annot_df,
                       annotation_col = gene_annot_df,
                       annotation_colors = annot_col)
g <- pheatmap2edit(g)
g <-pheatmap::pheatmap(check, color = my_color,cluster_rows = F,cluster_cols = F,
                       breaks = my_breaks,border_color = NA,annotation_row = gene_annot_df,
                       na_col = "white",
                       annotation_col = gene_annot_df,
                       annotation_colors = annot_col)

g <- pheatmap2edit(g)
g