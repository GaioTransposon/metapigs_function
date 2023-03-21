library(readr)
library(dplyr)
library(ggplot2)
library(tidyverse)

egg_all <- read_csv("/Users/dgaio/contigs/prodigal/eggnogg/KEGG/all_pathway_ko00500.csv") 
rec_all <- read_csv("/Users/dgaio/contigs/prodigal/reCOGnizer_results/KEGG/all_pathway_ko00500.csv") 

head(egg_all)
##
egg <- egg_all %>% group_by(pig,date,KO) %>% dplyr::summarise(sum=mean(norm_mapped_wa))
rec <- rec_all %>% group_by(pig,date,KO) %>% dplyr::summarise(sum=mean(norm_mapped_wa))

df <- merge(rec,egg, by=c("pig","date","KO"), all = TRUE)

df %>%
  ggplot(., aes(x=log(sum.x), y=log(sum.y))) +
  geom_point() +
  stat_smooth(method = "lm", col = "red") 
##



find_and_concat_fc <- function(path_to_fold_changes_files, fc) {
  
  these_files <- list.files(path = path_to_fold_changes_files, pattern=fc)
  mylist=list()
  for (i in seq_along(these_files)) {
    f_path=paste0(path_to_fold_changes_files,these_files[i])
    mylist[[i]] <- read.csv(file = f_path) %>%
      dplyr::filter(!significance=="ns")
  }
  # concantenate 
  only_significant <- do.call(rbind, mylist) 
  
  return(only_significant)
}

rec <- find_and_concat_fc("/Users/dgaio/contigs/prodigal/reCOGnizer_results/KEGG/","fc_t0_t10")
egg <- find_and_concat_fc("/Users/dgaio/contigs/prodigal/eggnogg/KEGG/", "fc_t0_t10")
merged <- merge(rec,egg,by=c("KO","pathway"), all=TRUE)
merged %>%
  ggplot(., aes(x=log_fc.x, y=log_fc.y)) +
  geom_point() +
  stat_smooth(method = "lm", col = "red") +
  labs(x="reCOGnizer annotations (log_fc)",
       y="eggnog annotations (log_fc)",
       title = "interval t0-t10")


rec <- find_and_concat_fc("/Users/dgaio/contigs/prodigal/reCOGnizer_results/KEGG/","fc_t2_t8")
egg <- find_and_concat_fc("/Users/dgaio/contigs/prodigal/eggnogg/KEGG/", "fc_t2_t8")
merged <- merge(rec,egg,by=c("KO","pathway"), all=TRUE)
merged %>%
  ggplot(., aes(x=log_fc.x, y=log_fc.y)) +
  geom_point() +
  stat_smooth(method = "lm", col = "red") +
  labs(x="reCOGnizer annotations (log_fc)",
       y="eggnog annotations (log_fc)",
       title = "interval t2-t8")




shift_agreement_check <- function(path_to_fold_changes_files) {
  
  t0t10 <- list.files(path = path_to_fold_changes_files, pattern="fc_t0_t10")
  mylist=list()
  for (i in seq_along(t0t10)) {
    f_path=paste0(path_to_fold_changes_files,t0t10[i])
    mylist[[i]] <- read.csv(file = f_path) %>%
      dplyr::filter(!significance=="ns")
  }
  # concantenate 
  t0t10_sign <- do.call(rbind, mylist) 
  
  t2t8 <- list.files(path = path_to_fold_changes_files, pattern="fc_t2_t8")
  mylist=list()
  for (i in seq_along(t2t8)) {
    f_path=paste0(my_path,t2t8[i])
    mylist[[i]] <- read.csv(file = f_path) %>%
      dplyr::filter(!significance=="ns") 
  }
  # concantenate 
  t2t8_sign <- do.call(rbind, mylist) 
  
  merged <- merge(t0t10_sign,t2t8_sign,by=c("KO","pathway"), all=TRUE)
  return(merged)
}


x <- shift_agreement_check("/Users/dgaio/contigs/prodigal/reCOGnizer_results/KEGG/")
x %>%
  ggplot(., aes(x=log_fc.x, y=log_fc.y)) +
  geom_point() +
  stat_smooth(method = "lm", col = "red") +
  labs(x="log fc - interval t0-t10",
       y="log fc - interval t0-t10",
       title = "reCOGnizer")


y <- shift_agreement_check("/Users/dgaio/contigs/prodigal/eggnogg/KEGG/")
merged %>%
  ggplot(., aes(x=log_fc.x, y=log_fc.y)) +
  geom_point() +
  stat_smooth(method = "lm", col = "red") +
  labs(x="log fc - interval t0-t10",
       y="log fc - interval t0-t10",
       title = "eggnog")









path_to_pathways="/Users/dgaio/contigs/prodigal/eggnogg/KEGG/"
these_files <- list.files(path = path_to_pathways, pattern="^all_")

mylist=list()
for (i in seq_along(these_files[1:3])) {
  f_path=paste0(path_to_pathways,these_files[i])
  keep <- read.csv(file = f_path) %>%
    # first keep only data of subjects that have been sampled at all time points:  
    dplyr::filter(date=="t0"|date=="t2"|date=="t4"|date=="t6"|date=="t8"|date=="t10") %>%
    dplyr::select(pig,date) %>% distinct() %>% group_by(pig) %>% tally() %>% dplyr::filter(n==6) %>% dplyr::select(pig)

}







egg_all <- read_csv("/Users/dgaio/contigs/prodigal/eggnogg/KEGG/all_pathway_ko00500.csv") 

View(egg_all)
# first keep only data of subjects that have been sampled at all time points:  
keep <- egg_all %>%
  dplyr::filter(date=="t0"|date=="t2"|date=="t4"|date=="t6"|date=="t8"|date=="t10") %>%
  dplyr::select(pig,date) %>% distinct() %>% group_by(pig) %>% tally() %>% dplyr::filter(n==6) %>% dplyr::select(pig)
egg_sub <- egg_all[egg_all$pig %in% keep$pig,]
# filter main time points: t0, t2, t4, t6, t8, t10. 
df <- egg_sub %>%
  dplyr::filter(date=="t0"|date=="t2"|date=="t4"|date=="t6"|date=="t8"|date=="t10")
# order:
df$date <- factor(df$date,levels=c("t0","t2","t4","t6","t8","t10"))

df <- df %>%
  group_by(pig,KO,date) %>%
  dplyr::summarise(x=mean(norm_mapped_wa), .groups = "drop") 


# split 
multiple_DFs <- split( df , f = df$KO ,drop = TRUE)
NROW(multiple_DFs)









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



