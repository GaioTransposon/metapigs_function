
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(patchwork)
# downloaded ComplexHeatmamp locally from: https://cran.r-project.org/src/contrib/Archive/rjson/
library(dendextend)
library(ComplexHeatmap)
library("cluster")
library(scales)

dfs <- read.csv(file = '/Users/dgaio/contigs/prodigal/eggnogg/KEGG/heatmap.csv')

# dfs <- dfs %>%
#   dplyr::select(pathway_description,pathway,KO,t2,t8)

first <- as.list(dfs %>% 
                   group_by(pathway) %>%
                   tally() %>%
                   arrange(desc(n)) %>%
                   top_n(n = 4))

dfs_sub <- dfs[dfs$pathway %in% first$pathway,]

rownames(dfs_sub) <- paste0(dfs_sub$KO,'_',dfs_sub$pathway)
mylabels <- dfs_sub$pathway_description



dfs_sub$pathway_description <- NULL
dfs_sub$pathway <- NULL
dfs_sub$KO <- NULL


m <- as.matrix(dfs_sub)
m <- rescale(m, to = c(-1, 1))    # doesn't change how results look; only the legend is clearer

# split by a vector specifying rowgroups
Heatmap(m, name = "avg_ab",
        split = mylabels, 
        row_names_gp = gpar(fontsize = 4), 
        cluster_row_slices = TRUE, 
        clustering_distance_rows = "euclidean",
        show_row_dend = FALSE,
        cluster_columns = FALSE, 
        width = unit(6, "cm"), 
        row_title_rot = 0, 
        column_names_rot = 0, gap = unit(0.05, "cm"),
        border = "black",
        row_title_gp = gpar(fontsize = 7))


wo_pseudoc
