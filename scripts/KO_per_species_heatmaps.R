
if (getwd()!="/Users/dgaio") {
  .libPaths( c( .libPaths(), "/shared/homes/152324/R/library") )
  .libPaths( c( .libPaths(), "/shared/homes/152324/miniconda3/envs/recognizer_env/lib/R/library") )
} else {
  message("Running on local")
}



library(readr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(stringr)
library(tidyr)
library(data.table)
library(tidyverse)
library(pheatmap)
library(ComplexHeatmap)
library(matrixStats)
# BiocManager::install("InteractiveComplexHeatmap")
library(InteractiveComplexHeatmap)
library(shiny)


# VISUALIZATION:


filename <- paste0(my_path,"KO_species/","KO_to_species_bestMAGs.csv")
all_dfs <- read_csv(filename, col_types = cols(species = col_character()))

pathways_lengths <- read_csv("github/metapigs_function/middle_dir/pathways_lengths.csv")


NROW(all_dfs)
NROW(pathways_lengths)
merged <- merge(all_dfs,pathways_lengths, by.x="paths",by.y="pathway_description")
NROW(merged)
head(merged)
unique(merged$paths)

merged <- merged %>%
  dplyr::filter(!map_ko=="ko00363") %>% # removing as this pathway is only made of 2 KOs (the next bigger one is made out of 7)
  dplyr::mutate(ratio=n/KOs_per_pathway) %>%
  dplyr::select(paths,species,ratio)
head(merged)


dfs_sub <- as.data.frame(merged %>%
                           #subset(., (paths %in% carb) | (paths %in% ene) | (paths %in% bile) | (paths %in% aa)) %>% # 
                           pivot_wider(names_from = species, values_from = ratio, values_fill = 0)) # fills with 0 if NA
head(dfs_sub)

# get rid of pathways not informative - low variance per row
dfs_sub$row_var = rowVars(as.matrix(dfs_sub[,-1]))
dfs_sub <- dfs_sub %>%
  dplyr::arrange(desc(row_var)) %>% # highest variance on top
  slice(1:ceiling(NROW(dfs_sub)/100*90)) # remove paths with lowest variance (10% of tot paths)

rownames(dfs_sub) <- dfs_sub$paths
mylabels <- dfs_sub$paths
dfs_sub$paths <- NULL
dfs_sub$row_var <- NULL



m <- as.matrix(dfs_sub)



# split by a vector specifying rowgroups
my_m <- Heatmap(m, name = "avg_ab",
                #row_names_gp = gpar(fontsize = 4),
                #column_names_gp = gpar(fontsize = 4),
                clustering_distance_rows = "euclidean",
                show_row_dend = FALSE,
                cluster_rows = TRUE,
                cluster_columns = TRUE,
                width = unit(20, "cm"),
                row_title_rot = 0,
                gap = unit(0.02, "cm"),
                border = "black",row_names_rot = 0,
                row_title_gp = gpar(fontsize = 7))


pheatmap(m, scale="column", cluster_rows = TRUE, cluster_cols= TRUE,
         fontsize_col = 1, fontsize_row = 4, treeheight_row = 1,
         color=colorRampPalette(c("navy", "white", "red"))(20))



ht1 = draw(my_m)

ui = fluidPage(
  h3("My first interactive ComplexHeatmap Shiny app"),
  p("This is an interactive heatmap visualization on a random matrix."),
  InteractiveComplexHeatmapOutput()
)

server = function(input, output, session) {
  makeInteractiveComplexHeatmap(input, output, session, ht1)
}

shinyApp(ui, server)



# #############
# m.pca2 <- prcomp(m[,-1], center = TRUE,scale. = TRUE)
# paths <- substr(m$paths, start = 1, stop = 8)
# 
# fviz_pca_ind(m.pca2,
#              geom.ind="point",
#              #fill.ind = dates, #col.ind = rainbow(n = 11),
#              pointshape = 21, pointsize = 2,
#              habillage = paths,
#              #geom.ind = "point", # show points only (nbut not "text")
#              col.ind = paths, # color by groups
#              #palette = c("#00AFBB", "#E7B800", "#FC4E07"),
#              addEllipses = FALSE, # Concentration ellipses
#              title="")+
#   scale_color_manual(name="time point",
#                      values=rainbow(n = NROW(m)))+
#   theme(legend.position="none")+
#   guides(color = guide_legend(nrow = 1))
# 
# fviz_pca_biplot(m.pca2,
#                 axes=c(1,2),   # which PCA ordinates
#                 select.var = list(cos2 = 20), # how many species to show - top n
#                 #select.ind = list(cos2 = 20), # how many paths to show -  top n
#                 label = "var",
#                 geom.ind=c("point","text"),
#                 geom.var=c("arrow","text"),
#                 pointshape = 21, pointsize = 2,
#                 habillage = paths,
#                 alpha.var ="contrib",
#                 col.ind = "#00AFBB",
#                 repel = TRUE,
#                 labelsize=1) +
#   scale_color_manual(name="paths",
#                      values=rainbow(n = NROW(m)))+
#   ggtitle("")+
#   theme(legend.position="right",
#         panel.border = element_rect(colour = "black", fill=NA, size=1))
# #############
