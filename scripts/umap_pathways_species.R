
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



# paths to use: 
my_path <- "~/contigs/prodigal/eggnogg/KEGG/"
print(my_path)

these_files <- list.files(path = my_path, pattern="^all")

these_files <- these_files

df_list=list()
for (f in these_files) {
  
  f_path=paste0(my_path,f)
  df <- read_csv(f_path, col_types = cols(species = col_character()))
  
  if (NROW(df) > 0) {
    
    print(f_path)
    
    descr <- unique(df$pathway_description)
    
    essential <- df %>%
      dplyr::select(KO,species) %>%
      distinct()
    
    # if species is NA, assign no_bin
    essential$species <- essential$species %>% replace_na('no_bin')
    
    save <- essential %>%
      group_by(species) %>%
      tally()
    
    df_list[[descr]] <- save
  
  }
    
}



all_dfs <- rbindlist(df_list, idcol = "paths") 


filename <- paste0(my_path,"KO_species/","KO_to_species.csv")
fwrite(all_dfs, file=filename)
head(all_dfs)
NROW(all_dfs)









# 
# library(scales)
# library(magrittr)
# library(pacman)
# library(textshape)
# 
# 
# 
# # carb,ene,lpid,bile,nuc,aa,oth_aa,glycan
# #all_dfs <- subset(all_dfs, (paths %in% carb) | (paths %in% ene) | (paths %in% bile) | (paths %in% aa))
# 
# 
#                      
#                      
#                      
# m <- all_dfs %>%                     
#   pivot_wider(names_from = species, values_from = n, values_fill = 0) # fills with 0 if NA
# head(m)
# class(m)
# 
# # transform to percentages
# m[,-1] <- lapply(m[,-1], function(x) round(x / sum(x)*100,digits = 2))
# head(m)
# 
# 
# # rownames(m) <- m$paths
# # rownames(m)
# # m$paths <- NULL
# 
# rownames(m)
# 
# head(m)
# 
# m <- cluster_matrix(m, dim = "row")
# 
# head(m)
# 
# 
# 
# 
# library(tsne)
# library(plotly)
# 
# # iris$Species == m$paths
# 
# features <- subset(m[,-1]) 
# 
# set.seed(0)
# tsne <- tsne(features, initial_dims = 3, k =3)
# tsne <- data.frame(tsne)
# pdb <- cbind(tsne,m$paths)
# options(warn = -1)
# fig <-  plot_ly(data = pdb ,x =  ~X1, y = ~X2, z = ~X3, color = ~m$paths, ) %>% 
#   layout()
#   # add_markers(size = 8) %>%
#   # layout( 
#   #   xaxis = list(
#   #     zerolinecolor = "#ffff",
#   #     zerolinewidth = 2,
#   #     gridcolor='#ffff'), 
#   #   yaxis = list(
#   #     zerolinecolor = "#ffff",
#   #     zerolinewidth = 2,
#   #     gridcolor='#ffff'),
#   #   scene =list(bgcolor = "#e5ecf6"))
# fig

