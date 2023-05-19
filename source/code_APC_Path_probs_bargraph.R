.libPaths(Sys.getenv("Yoyo_revised"))
sapply(list.files("D:/OneDrive - dgist.ac.kr/Function/",full.names = T), source)
library(Seurat)
library(dplyr)
library(ggplot2)

APC_seurat <- readRDS("D:/OneDrive - dgist.ac.kr/ADIPO_crtc/KU_CRTC/data/seurat/qc_positive_APC_seurat.rds")

probs_df_stage2 <- data.frame()
stage2_seurat <- APC_seurat[,APC_seurat$apc_annot_simple=="APC.stage2"]
for( i in levels(stage2_seurat$condition)){
  path1_probs <-c(i,"Path1",mean(stage2_seurat$Path1[stage2_seurat$condition==i]),
                  sd(stage2_seurat$Path1[stage2_seurat$condition==i]), 
                  ncol(stage2_seurat[,stage2_seurat$condition==i]))
  path2_probs <-c(i,"Path2",mean(stage2_seurat$Path2[stage2_seurat$condition==i]),
                  sd(stage2_seurat$Path2[stage2_seurat$condition==i]),
                  ncol(stage2_seurat[,stage2_seurat$condition==i]))
  tmp_df <- rbind(path1_probs,path2_probs)
  probs_df_stage2 <- rbind(probs_df_stage2,tmp_df)
}

colnames(probs_df_stage2) <- c("Condition","Path","Mean","SD","n")

df <- data.frame(row.names = 1:(1*ncol(stage2_seurat)))
df$probs <- c(stage2_seurat$Path1)
df$path <- c(rep("path1",ncol(stage2_seurat)))
df$condition <- as.character(stage2_seurat$condition)

# Box plots with jittered points
library(ggpubr)
p1 <- ggboxplot(df, x = "condition", y = "probs",
                color = "condition", palette = c("#4374D9","#3DB7CC","#CC3D3D","#CC723D"),
                 shape = "condition")
# Add p-values comparing groups
my_comparisons <- list( c("Young_WT", "Young_KO"), c("Young_WT", "Old_WT"), c("Old_WT", "Old_KO"))
p2 <- p1 + stat_compare_means(comparisons = my_comparisons,
                              symnum.args=list(cutpoints = c(0, 0.01, 0.05, 1), symbols = c( "**", "*", "ns")))


