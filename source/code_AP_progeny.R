library(progeny)
library(dplyr)
library(Seurat)
library(ggplot2)
library(tidyr)
library(readr)
library(pheatmap)
library(tibble)
APC_seurat <- readRDS("qc_positive_APC_seurat.rds")
Idents(APC_seurat) <- factor(paste0(APC_seurat$apc_annot_simple,"_",APC_seurat$condition),
                             levels = paste0(rep(c("AP1","AP2"),each=4)
                                             ,"_",rep(levels(APC_seurat$condition),2)))
# APC_seurat <- subset(APC_seurat,idents = c(1,2))

APC_seurat <- ScaleData(APC_seurat, features = rownames(APC_seurat))
CellsClusters <- data.frame(Cell = names(Idents(APC_seurat)), 
                            CellType = Idents(APC_seurat),
                            stringsAsFactors = FALSE)



## We compute the Progeny activity scores and add them to our Seurat object
## as a new assay called Progeny.
APC_seurat <- progeny(APC_seurat, scale=FALSE, organism="Mouse",top=500, perm=1,
                    return_assay = TRUE)
## We can now directly apply Seurat functions in our Progeny scores. 
## For instance, we scale the pathway activity scores. 
APC_seurat <- Seurat::ScaleData(APC_seurat, assay = "progeny") 
## We transform Progeny scores into a data frame to better handling the results
progeny_scores_df <- 
  as.data.frame(t(GetAssayData(APC_seurat, slot = "scale.data", 
                               assay = "progeny"))) %>%
  rownames_to_column("Cell") %>%
  gather(Pathway, Activity, -Cell) 
## We match Progeny scores with the cell clusters.
progeny_scores_df <- inner_join(progeny_scores_df, CellsClusters)

## We summarize the Progeny scores by cellpopulation
summarized_progeny_scores <- progeny_scores_df %>% 
  group_by(Pathway, CellType) %>%
  summarise(avg = mean(Activity), std = sd(Activity))

## We prepare the data for the plot
summarized_progeny_scores_df <- summarized_progeny_scores %>%
  dplyr::select(-std) %>%   
  spread(Pathway, avg) %>%
  data.frame(row.names = 1,  stringsAsFactors = FALSE) 

paletteLength = 100
library(RColorBrewer)
paletteLength = 100
myColor = colorRampPalette(c("Darkblue", "white","red"))(paletteLength)
#my_color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(paletteLength)

progenyBreaks = c(seq(min(summarized_progeny_scores_df), 0, 
                      length.out=ceiling(paletteLength/2) + 1),
                  seq(max(summarized_progeny_scores_df)/paletteLength, 
                      max(summarized_progeny_scores_df), 
                      length.out=floor(paletteLength/2)))
source("D:/OneDrive - dgist.ac.kr/Function/create_pptx_function.R")
progeny_hmap = pheatmap(t(summarized_progeny_scores_df),fontsize=14, 
                        cluster_cols = F,
                        fontsize_row = 10, 
                        color=myColor, breaks = progenyBreaks, 
                        main = "PROGENy (500)", angle_col = 45,
                        border_color = NA, treeheight_row = 0)
progeny_hmap <- pheatmap2edit(progeny_hmap)
progeny_hmap = pheatmap(t(summarized_progeny_scores_df[,c("Hypoxia","TGFb","WNT","PI3K","VEGF","MAPK")]),fontsize=14, 
                        cluster_cols = F,
                        fontsize_row = 10, 
                        color=myColor, breaks = progenyBreaks, 
                        main = "PROGENy (500)", angle_col = 45,
                        border_color = NA, treeheight_row = 0)
progeny_hmap <- pheatmap2edit(progeny_hmap)
