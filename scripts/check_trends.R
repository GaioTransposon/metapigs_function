library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(patchwork)
library("cluster")
library(scales)
library(stringr)
library(readr)
library(ggpubr)
library(EnvStats)
library(purrr)
library(data.table)
library(tidyverse)


ko00071 <- read.csv(file = '/Users/dgaio/contigs/prodigal/eggnogg/KEGG/all_pathway_ko00071.csv')
NROW(ko00071)

# Filter out 1 ortholog gene that shows a time trend, for example: K00022
K00022 <- ko00071 %>%
  dplyr::filter(KO=="K00022") 

K00022["species"][K00022["species"] == ''] <- "none"
K00022$species

# reorder dates 
K00022$date  = factor(K00022$date, levels=c("t0",
                                            "t1", 
                                            "t2",
                                            "t3",
                                            "t4",
                                            "t5",
                                            "t6",
                                            "t7",
                                            "t8",
                                            "t9",
                                            "t10",
                                            "tM"))


# subset to subjects that have been sampled till the end 
sampled <- as.list(K00022 %>%
                     dplyr::filter(date=="t9") %>%
                     dplyr::select(pig) %>%
                     distinct())
NROW(K00022)
K00022 <- subset(K00022, (pig %in% sampled$pig))
NROW(K00022)



# Prevalence among subjects? 
n_pigs <- K00022 %>%
  dplyr::select(pig) %>%
  distinct() %>%
  tally()

# Which species bear it? 
n_species <- K00022 %>%
  dplyr::select(species) %>%
  distinct() %>%
  tally()


pseudo <- min(K00022$norm_mapped_wa[K00022$norm_mapped_wa > 0])
p1 <- K00022 %>%
  dplyr::mutate(norm_mapped_wa=norm_mapped_wa+pseudo) %>%
  group_by(pig,date) %>%
  dplyr::summarise(norm_mapped_wa=mean(norm_mapped_wa)) %>%
  ggplot(., aes(x=date,y=log(norm_mapped_wa))) +
  geom_boxplot()+
  labs(#title = KO,
       #subtitle = paste0(KO,"_",funct_cat),
       caption = as.character(paste0("tot# subjects:",n_pigs,
                                     "\ntot# species:",n_species-1))) + # min 1 because otherwise NA is counted
  stat_n_text(hjust=0, size = 2, angle = 90)+
  theme(title = element_text(size=5))


# top 10 species and plot 
these_species <- K00022 %>%
  group_by(species) %>%
  tally() %>%
  dplyr::mutate(perc=round(n/sum(n)*100,2)) %>%
  dplyr::arrange(desc(perc)) %>%
  head(10)

df_temp <- K00022 %>%
  #dplyr::mutate(norm_mapped_wa=norm_mapped_wa+pseudo) %>%
  group_by(pig,species,date) %>%
  dplyr::summarise(norm_mapped_wa=mean(norm_mapped_wa)) 


plot_species <- left_join(these_species, df_temp)
plot_species$species_perc=paste0(plot_species$species,"\n",plot_species$perc)


p2 <- plot_species %>%
  ggplot(., aes(x=date,y=log(norm_mapped_wa))) +
  geom_boxplot()+
  labs(#title = KO,
    #subtitle = paste0(KO,"_",funct_cat),
    caption = as.character(paste0("tot# subjects:",n_pigs,
                                  "\ntot# species:",n_species-1))) + # min 1 because otherwise NA is counted
  theme(title = element_text(size=5))+
  facet_wrap(~factor(species_perc, levels=unique(plot_species$species_perc)))+
  labs(title = plot_species$perc[0]) +
  stat_n_text(hjust=0, size = 1, angle = 90)+
  theme(strip.text = element_text(size=5))


both <- ggarrange(
  p1,p2,ncol=2)


          
         
     















