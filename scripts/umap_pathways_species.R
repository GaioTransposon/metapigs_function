
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










library(scales)
library(magrittr)
library(pacman)
library(textshape)




m <- all_dfs %>%
  subset(., (paths %in% carb) | (paths %in% ene) | (paths %in% bile) | (paths %in% aa) ) %>% # carb,ene,lpid,bile,nuc,aa,oth_aa,glycan
  pivot_wider(names_from = species, values_from = n, values_fill = 0) # fills with 0 if NA
head(m)
class(m)


# transform to percentages
m[,-1] <- lapply(m[,-1], function(x) round(x / sum(x)*100,digits = 2))
head(m)


# rownames(m) <- m$paths
# rownames(m)
# m$paths <- NULL

rownames(m)

head(m)

m <- cluster_matrix(m, dim = "row")

head(m)




NROW(m)





library(palmerpenguins)
library(umap)
theme_set(theme_bw(18))

penguins <- penguins %>% 
  drop_na() %>%
  select(-year)%>%
  mutate(ID=row_number()) 
penguins_meta <- penguins %>%
  select(ID, species, island, sex)

set.seed(142)
umap_fit <- penguins %>%
  select(where(is.numeric)) %>%
  column_to_rownames("ID") %>%
  scale() %>% 
  umap()

umap_df <- umap_fit$layout %>%
  as.data.frame()%>%
  rename(UMAP1="V1",
         UMAP2="V2") %>%
  mutate(ID=row_number())%>%
  inner_join(penguins_meta, by="ID")

umap_df %>% head()

umap_df %>%
  ggplot(aes(x = UMAP1, 
             y = UMAP2, 
             color = species,
             shape = sex))+
  geom_point()+
  labs(x = "UMAP1",
       y = "UMAP2",
       subtitle = "UMAP plot")

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

head(m)


# m[,1:30] %>% pivot_longer(cols=-paths) %>%
#   ggplot(., aes(x=reorder(name,value, mean),y=value,group=paths, fill=paths))+
#   geom_bar(stat = "identity") 



library(tsne)
library(plotly)

# iris$Species == m$paths

features <- subset(m[,-1])

set.seed(0)
tsne <- tsne(features, initial_dims = 3, k =3)
tsne <- data.frame(tsne)
pdb <- cbind(tsne,m$paths)
options(warn = -1)
fig <-  plot_ly(data = pdb ,x =  ~X1, y = ~X2, z = ~X3, color = ~m$paths) %>%
  add_markers(size = 8) %>%
  layout(
    xaxis = list(
      zerolinecolor = "#ffff",
      zerolinewidth = 2,
      gridcolor='#ffff'),
    yaxis = list(
      zerolinecolor = "#ffff",
      zerolinewidth = 2,
      gridcolor='#ffff'),
    scene =list(bgcolor = "#e5ecf6"))
fig

