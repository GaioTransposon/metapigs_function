
# #############
# # n<-t(m)
# # rownames(n)
# m.pca2 <- prcomp(m, center = TRUE,scale. = TRUE)
# paths <- substr(rownames(m), start = 1, stop = 8)
# 
# fviz_pca_ind(m.pca2,
#              #fill.ind = dates, #col.ind = rainbow(n = 11),
#              pointshape = 21, pointsize = 2,
#              habillage = paths,
#              geom.ind = "point", # show points only (nbut not "text")
#              col.ind = paths, # color by groups
#              #palette = c("#00AFBB", "#E7B800", "#FC4E07"),
#              addEllipses = FALSE, # Concentration ellipses
#              title="")+
#   scale_color_manual(name="time point",
#                      values=rainbow(n = NROW(m)))+
#   theme(legend.position="bottom")+
#   guides(color = guide_legend(nrow = 1))
# 
# fviz_pca_biplot(m.pca2,
#                 axes=c(1,2),   # which PCA ordinates
#                 #select.var = list(cos2 = 20), # how many species to show - top n
#                 select.ind = list(cos2 = 10), # how many paths to show -  top n
#                 #label = "var",
#                 #geom.ind=c("point","text"),
#                 geom.var=c("point"),
#                 pointshape = 21, pointsize = 2,
#                 habillage = paths,
#                 #alpha.var ="contrib",
#                 #col.ind = "#00AFBB",
#                 repel = TRUE,
#                 labelsize=1) +
#   scale_color_manual(name="paths",
#                      values=rainbow(n = NROW(m)))+
#   ggtitle("")+
#   theme(legend.position="right",
#         panel.border = element_rect(colour = "black", fill=NA, size=1))
# #############


# 
# library(tsne)
# library(plotly)
# data("iris")
# 
# 
# dfs_sub <- as.data.frame(merged %>%
#                            #subset(., (paths %in% carb) | (paths %in% ene) | (paths %in% bile) | (paths %in% aa)) %>% # 
#                            pivot_wider(names_from = species, values_from = ratio, values_fill = 0)) # fills with 0 if NA
# head(dfs_sub)
# 
# # get rid of pathways not informative - low variance per row
# dfs_sub$row_var = rowVars(as.matrix(dfs_sub[,-1]))
# dfs_sub <- dfs_sub %>%
#   dplyr::arrange(desc(row_var)) %>% # highest variance on top
#   slice(1:ceiling(NROW(dfs_sub)/100*90)) %>% # remove paths with lowest variance (10% of tot paths)
#   dplyr::select(!row_var)
# 
# features <- subset(dfs_sub, select = -c(paths)) 
# 
# set.seed(0)
# tsne <- tsne(features, initial_dims = 2)
# tsne <- data.frame(tsne)
# pdb <- cbind(tsne,dfs_sub$paths)
# options(warn = -1)
# fig <-  plot_ly(data = pdb ,x =  ~X1, y = ~X2, 
#                 type = 'scatter', 
#                 mode = 'markers', split = ~dfs_sub$paths)
# 
# fig <- fig %>%
#   layout(
#     plot_bgcolor = "#e5ecf6"
#   )
# 
# fig
# 
# 
# #set.seed(0)
# features <- subset(dfs_sub, select = -c(paths)) 
# tsne <- tsne(features, initial_dims = 3, k =3)
# tsne <- data.frame(tsne)
# pdb <- cbind(tsne,dfs_sub$paths)
# options(warn = -1)
# fig <-  plot_ly(data = pdb ,x =  ~X1, y = ~X2, z = ~X3, color = ~dfs_sub$paths) %>% 
#   add_markers(size = 8) %>%
#   layout( 
#     xaxis = list(
#       zerolinecolor = "#ffff",
#       zerolinewidth = 2,
#       gridcolor='#ffff'), 
#     yaxis = list(
#       zerolinecolor = "#ffff",
#       zerolinewidth = 2,
#       gridcolor='#ffff'),
#     scene =list(bgcolor = "#e5ecf6"))
# fig




# 
# 
# library(plotly) 
# library(umap) 
# dfs_sub.data = dfs_sub[,-1] 
# dfs_sub.labels = dfs_sub[, "paths"] 
# dfs_sub.umap = umap(dfs_sub.data, n_components = 2, random_state = 15) 
# layout <- dfs_sub.umap[["layout"]] 
# layout <- data.frame(layout) 
# final <- cbind(layout, dfs_sub$paths) 
# 
# fig <- plot_ly(final, x = ~X1, y = ~X2, color = ~dfs_sub$paths, type = 'scatter', mode = 'markers')%>%  
#   layout(
#     plot_bgcolor = "#e5ecf6",
#     legend=list(title=list(text='species')), 
#     xaxis = list( 
#       title = "0"),  
#     yaxis = list( 
#       title = "1")) 
# 
# dfs_sub.umap = umap(dfs_sub.data, n_components = 3, random_state = 15) 
# layout <- dfs_sub.umap[["layout"]] 
# layout <- data.frame(layout) 
# final <- cbind(layout, dfs_sub$paths) 
# 
# fig2 <- plot_ly(final, x = ~X1, y = ~X2, z = ~X3, color = ~dfs_sub$paths) 
# fig2 <- fig2 %>% add_markers() 
# fig2 <- fig2 %>% layout(scene = list(xaxis = list(title = '0'), 
#                                      yaxis = list(title = '1'), 
#                                      zaxis = list(title = '2'))) 
# 
# fig2






library(plotly) 
library(umap) 




dfs_sub <- as.data.frame(merged %>%
                           #subset(., (paths %in% carb) | (paths %in% ene) | (paths %in% bile) | (paths %in% aa)) %>% # 
                           pivot_wider(names_from = paths, values_from = ratio, 
                                       values_fill = 0, names_prefix = "path_")) # fills with 0 if NA
head(dfs_sub)

# get rid of non informative species not informative - low variance per row
dfs_sub$row_var = rowVars(as.matrix(dfs_sub[,-1]))
# dfs_sub <- dfs_sub %>%
#   dplyr::arrange(desc(row_var)) %>% # highest variance on top
#   slice(1:ceiling(NROW(dfs_sub)/100*90)) %>% # remove species with lowest variance (10% of tot species)
#   dplyr::select(!row_var)
NROW(dfs_sub)
dfs_sub <- inner_join(dfs_sub,gtdb)
NROW(dfs_sub)


head(dfs_sub)

dfs_sub.data = dfs_sub %>% select(starts_with("path"))

dfs_sub.labels = dfs_sub[, "species"] 
dfs_sub.umap = umap(dfs_sub.data, n_components = 3, random_state = 15) 
layout <- dfs_sub.umap[["layout"]] 
layout <- data.frame(layout) 
final <- cbind(layout, dfs_sub$species) 

plot_ly(final, x = ~X1, y = ~X2, z = ~X3,
        split = ~dfs_sub$class,
        legendgroup = ~dfs_sub$species, 
        name = ~dfs_sub$class, 
        color = ~dfs_sub$class,
        hovertext = ~ paste0(dfs_sub$family, "__", dfs_sub$species),
        hoverinfo = "text",
        type = 'scatter', mode = 'markers') %>%  
  layout(
    plot_bgcolor = "#e5ecf6",
    legend=list(title=list(text='Class')),
    xaxis = list(
      title = "0"),
    yaxis = list(
      title = "1"))


# dfs_sub.umap = umap(dfs_sub.data, n_components = 3, random_state = 15) 
# layout <- dfs_sub.umap[["layout"]] 
# layout <- data.frame(layout) 
# final <- cbind(layout, dfs_sub$species) 
# 
# fig2 <- plot_ly(final, x = ~X1, y = ~X2, z = ~X3, color = ~dfs_sub$family) 
# fig2 <- fig2 %>% add_markers() 
# fig2 <- fig2 %>% layout(scene = list(xaxis = list(title = '0'), 
#                                      yaxis = list(title = '1'), 
#                                      zaxis = list(title = '2'))) 





# gtdb inventing names? 



gtdb_bins_completeTaxa <- read_csv("~/github/metapigs_dry/middle_dir/gtdb_bins_completeTaxa", 
                                   col_types = cols(pig = col_character()))
gtdb <- as.data.frame(gtdb_bins_completeTaxa %>% 
                        dplyr::select(domain,phylum,class,order,family,genus,species) %>%
                        distinct())



library(umap)



dfs_sub <- as.data.frame(merged %>%
                           #subset(., (paths %in% carb) | (paths %in% ene) | (paths %in% bile) | (paths %in% aa)) %>% # 
                           pivot_wider(names_from = paths, values_from = ratio, 
                                       values_fill = 0, names_prefix = "path_")) # fills with 0 if NA
head(dfs_sub)

# get rid of non informative species not informative - low variance per row
dfs_sub$row_var = rowVars(as.matrix(dfs_sub[,-1]))
# dfs_sub <- dfs_sub %>%
#   dplyr::arrange(desc(row_var)) %>% # highest variance on top
#   slice(1:ceiling(NROW(dfs_sub)/100*90)) %>% # remove species with lowest variance (10% of tot species)
#   dplyr::select(!row_var)
NROW(dfs_sub)
dfs_sub <- inner_join(dfs_sub,gtdb)
NROW(dfs_sub)


head(dfs_sub)

dfs_sub.data = dfs_sub %>% select(starts_with("path"))

dfs_sub.labels = dfs_sub[, "species"] 
dfs_sub.umap = umap(dfs_sub.data, n_components = 2, random_state = 15) 
layout <- dfs_sub.umap[["layout"]] 
layout <- data.frame(layout) 
final <- cbind(layout, dfs_sub$species) 

plot_ly(final, x = ~X1, y = ~X2, 
        split = ~dfs_sub$phylum,
        legendgroup = ~dfs_sub$species, 
        name = ~dfs_sub$phylum, 
        color = ~dfs_sub$phylum,
        hovertext = ~ paste0(dfs_sub$family, "__", dfs_sub$species),
        hoverinfo = "text",
        type = 'scatter', mode = 'markers') %>%  
  layout(
    plot_bgcolor = "#e5ecf6",
    legend=list(title=list(text='Phylum')),
    xaxis = list(
      title = "0"),
    yaxis = list(
      title = "1"))

