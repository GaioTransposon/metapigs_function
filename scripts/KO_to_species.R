
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


# paths to use: 
my_path <- "~/contigs/prodigal/eggnogg/KEGG/"
print(my_path)

these_files <- list.files(path = my_path, pattern="^all")


KO_list=list()
for (f in these_files) {
  
  f_path=paste0(my_path,f)
  df <- read_csv(f_path, col_types = cols(species = col_character()))
  print(f_path)
  
  file_name <- paste0(my_path,'KO_species/',str_replace(f,".csv","_KO_species.pdf"))

  essential <- df %>%
    dplyr::select(pig,contig_orf,KO,species, pathway_description) %>%
    distinct()
  
  # if species is NA, assign no_bin
  essential["species"][essential["species"] == ''] <- "no_bin"
  essential$species <- essential$species %>% replace_na('no_bin')

  x <- essential %>%
    group_by(pathway_description,KO,species) %>%
    tally() %>%
    group_by(KO) %>%
    dplyr::mutate(perc=round(n/sum(n)*100,2)) %>%
    dplyr::arrange(KO,desc(perc)) %>%
    slice_head(n = 20)
  
  NROW(x)
  
  # save to file? 

  p1 <- x %>%
    ggplot(., aes(x=reorder(KO,perc, mean), y = perc, 
                  fill = species, group=species)) +
    geom_bar(stat = "identity", colour="black") +
    theme(legend.position = "no_bin",
          axis.text.x=element_text(size=5, angle=90))+
    ggtitle(label = x$pathway_description[1])
  
  # remove the KOs which are >=80% of the times found in unbinned contigs (not worth showing)
  remove <- x %>%
    dplyr::filter(perc>80  & species=="no_bin") %>%
    dplyr::select(KO)
  keep <- subset(x, !(KO %in% remove$KO))

  NROW(unique(keep$KO))
  # display for each, the top 3 species. 
  keep <- keep %>% 
    group_by(KO) %>%
    top_n(n = 3, wt = perc) %>%
    slice(1:3) 
  
  l=unique(keep$KO)
  chunk_length=30
  chunks <- split(l,ceiling(seq_along(l) / chunk_length))
  
  n=0
  plot_list = list()
  for (i in chunks) {
    n=n+1
    df_i <- subset(keep, (KO %in% i))
    plot_list[[n]] <- df_i %>%
      #drop.levels() %>% 
      ggplot(aes(x = reorder(KO, perc, FUN = mean), y = reorder(species, perc, FUN = mean), fill = perc)) +
      geom_tile() +
      scale_fill_distiller(type = "div", palette = "Spectral") +
      geom_text(aes(y=species,label=round(perc)), size=2)+
      facet_wrap(~KO, scales = "free",ncol = 3)+
      theme(axis.title.x=element_blank(),axis.text.x = element_blank(),
            axis.text.y=element_text(size=5),
            axis.title.y=element_blank(),
            strip.text.x = element_text(size = 7, colour = "black"),
            legend.key.height = unit(.5, "cm"),
            legend.key.width = unit(.5, "cm"),
            legend.key.size = unit(.5, "cm"),
            legend.text.align = 1,
            legend.title = element_text(size=7),
            legend.text = element_text(size=6),
            legend.position="bottom",
            legend.justification="bottom",
            legend.margin=margin(0,0,0,0),
            legend.box.margin=margin(-6,-6,-6,-6), # get closer farther from zero to get it closer to edge 
            axis.ticks.x = element_blank())+
      labs(fill="percentage of species carrying KO")+
      scale_y_discrete(labels=function(x){sub("\\s", "\n", x)})
  }
  
  pdf(file_name)
  print(p1)
  if (NROW(plot_list) !=0) {
    for (i in 1:NROW(plot_list)) {
      print(plot_list[[i]])
    }
  }
  dev.off()
  
}




