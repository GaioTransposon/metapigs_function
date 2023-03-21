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






