library(readr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(EnvStats)
library(purrr)
library(data.table)
library(tidyverse)


# paths to use: 
my_path <- paste0(getwd(),"/reCOGnizer_results/") # UTS HPC
# my_path="Desktop/contigs/prodigal/reCOGnizer_results/"    #local
print(my_path)

these_files <- list.files(path = my_path, pattern="significant.csv")

#these_files<-these_files[1:2]

for (f in these_files) {
  
  f_path=paste0(my_path,f)
  df <- read_csv(f_path)
  
  fname_new=str_replace(f, ".csv","_")
  fnew_path=paste0(my_path,fname_new)
  
  to_plot <- split( df , f = df$`CDD ID` )
  #NROW(to_plot)
  
  # splitting every 100 to print new pdf every 100 plots 
  seqq <- seq(from = 1, to = NROW(to_plot), by = 1)
  my_groups <- split(seqq, ceiling(seq_along(seqq)/100))
  
  counter=0
  for (i in my_groups) {
    
    counter=counter+1
    pdf_name=paste0(fnew_path,counter,'.pdf')
    
    # print new pdf every 100 plots 
    pdf(pdf_name)
    for (ii in i) {
      
      z <- to_plot[ii]
      z <- do.call(rbind.data.frame, z)
      funct_cat=as.character(z$`Functional category`[1])
      protein=as.character(z$`DB description`[1])
      CDD=as.character(z$`CDD ID`[1])
      
      n_pigs=z %>%
        dplyr::select(pig) %>%
        distinct() %>%
        tally()
      n_species=z %>%
        dplyr::select(species) %>%
        distinct() %>%
        tally()
      
      p1 <- z %>%
        ggplot(., aes(x=date,y=log(norm_mapped_wa))) +
        geom_boxplot()+
        labs(title = protein,
             subtitle = paste0(CDD,"_",funct_cat),
             caption = as.character(paste0("tot# subjects:",n_pigs,
                                           "\ntot# species:",n_species-1))) + # min 1 because otherwise NA is counted
        stat_n_text(vjust=-1)+
        theme(title = element_text(size=5))
      
      # top 10 species and plot 
      these_species <- z %>%
        group_by(species) %>%
        tally() %>%
        dplyr::mutate(perc=round(n/sum(n)*100,2)) %>%
        dplyr::arrange(desc(perc)) %>%
        head(10)
      
      z_sub <- inner_join(these_species,z)
      z_sub$species_perc=paste0(z_sub$species,"\n",z_sub$perc)
      
      # order the facets by number of species carrying the protein 
      p2 <- z_sub %>%
        ggplot(., aes(x=date,y=log(norm_mapped_wa))) +
        geom_boxplot()+
        facet_wrap(~factor(species_perc, levels=unique(z_sub$species_perc)))+
        labs(title = z_sub$perc[0]) +
        theme(strip.text = element_text(size=5))
      
      both <- ggarrange(
        p1,p2,ncol=2)
      
      plot(both)
      
    }
    
    dev.off()
  }
  
}













# #Pdf
# pdf('/Users/dgaio/Desktop/contigs/prodigal/reCOGnizer_results/test.pdf')
# #Loop
# for (i in 1:NROW(to_plot_sub)){
#   
#   z <- to_plot_sub[i]
#   z <- do.call(rbind.data.frame, z)
#   funct_cat=as.character(z$`Functional category`[1])
#   protein=as.character(z$`DB description`[1])
#   CDD=as.character(z$`CDD ID`[1])
#   
#   n_pigs=z %>%
#     dplyr::select(pig) %>%
#     distinct() %>%
#     tally()
#   n_species=z %>%
#     dplyr::select(species) %>%
#     distinct() %>%
#     tally()
# 
#   p1 <- z %>%
#     ggplot(., aes(x=date,y=log(norm_mapped_wa))) +
#     geom_boxplot()+
#     labs(title = name,
#          subtitle = paste0(CDD,"_",funct_cat),
#          caption = as.character(paste0("tot# subjects:",n_pigs,
#                                        "\ntot# species:",n_species-1))) + # min 1 because NA otherwise is counted
#     stat_n_text(vjust=-1)+
#     theme(title = element_text(size=5))
#   
#   # top 10 species and plot 
#   these_species <- z %>%
#     group_by(species) %>%
#     tally() %>%
#     dplyr::mutate(perc=round(n/sum(n)*100,2)) %>%
#     dplyr::arrange(desc(perc)) %>%
#     head(10)
#   
#   z_sub <- inner_join(these_species,z)
#   z_sub$species_perc=paste0(z_sub$species,"\n",z_sub$perc)
# 
#   # order the facets by number of species carrying the protein 
#   p2 <- z_sub %>%
#     ggplot(., aes(x=date,y=log(norm_mapped_wa))) +
#     geom_boxplot()+
#     facet_wrap(~factor(species_perc, levels=unique(z_sub$species_perc)))+
#     labs(title = z_sub$perc[0]) +
#     theme(strip.text = element_text(size=5))
#   
#   both <- ggarrange(
#     p1,p2,ncol=2)
# 
#   plot(both)
# 
# }
# dev.off()


















# 
# 
# pdf("/shared/homes/152324/contigs/prodigal/reCOGnizer_results/test2.pdf")
# recc %>%
#   dplyr::filter(h_statistic>500) %>%
#   #group_by(`CDD ID`) %>%
#   #tally()
#   group_by(pig, `CDD ID`, date) %>%
#   dplyr::summarise(tot = median(norm_mapped_wa, na.rm = TRUE)) %>%
#   group_by(`CDD ID`) %>%
#   dplyr::mutate(tot = tot/max(tot)) %>%
#   ungroup() %>%
#   ggplot(aes(x = factor(pig), y = reorder(`CDD ID`, tot, FUN = median), fill = tot)) +
#   geom_tile() +
#   scale_fill_distiller(type = "div", palette = "RdBu") +
#   facet_grid(~date, scales = "free") +
#   labs(x = "date", y = "CDD ID", fill = "normalized abundance")+
#   theme(axis.text.x=element_blank(),
#         axis.title.x=element_blank(),
#         legend.position="top")
# dev.off()
# 
# rec %>% 
#   dplyr::filter(h_statistic>100) %>%  #500
#   dplyr::select(`CDD ID`) %>%
#   distinct() %>%
#   tally()
#   
# 
# recc <- rec %>% 
#   dplyr::filter(h_statistic>100)
# 
# recc$pig_date=paste(recc$date, recc$pig, sep="_")
# 
# reccc <- recc %>%
#   group_by(pig_date,`CDD ID`) %>%
#   summarise(sum=sum(norm_mapped_wa)) %>%
#   pivot_wider(names_from = pig_date,values_from=sum) 
# 
# reccc <- as.data.frame(reccc)
# rownames(reccc) <- reccc$`CDD ID`  
# reccc$`CDD ID` <- NULL
# 
# ord=list("t0","t0","t0","t10","t10","t10")
# 
# 
# 
# library("gplots")
# heatmap.2(reccc, scale = "none", col = bluered(100), 
#           trace = "none", density.info = "none")
# 
# 
# 
# pdf("/shared/homes/152324/contigs/prodigal/reCOGnizer_results/test.pdf")
# pheatmap(reccc,clustering_distance_rows = "euclidean", cutree_rows = 8)  #scale="row" 
# pheatmap(reccc,clustering_distance_rows = "manhattan", scale="row", annotation_col = ord)  #scale="row" 
# pheatmap(reccc,clustering_distance_rows = "canberra")  #scale="row" 
# pheatmap(reccc,cluster_cols=FALSE,cluster_rows=TRUE)
# dev.off()
# 
# recc %>% 
#   dplyr::filter(`CDD ID`=="CDD:225343") %>%
#   ggplot(., aes(x=date,y=log(norm_mapped_wa))) +
#   geom_boxplot()+
#   stat_n_text(vjust=-1) +
#   facet_wrap(~pig) +
#   ggtitle(recc$`DB description`[1])
# 
# recc %>% 
#   dplyr::filter(`CDD ID`=="CDD:225343") %>%
#   group_by(species,pig) %>%
#   tally() %>%
#   arrange(desc(n))
# 
# View(recc)
# 
# ec <- recc %>%
#   dplyr::select(`EC number`) %>%
#   distinct()
# View(ec)
# getwd()
# write.table(ec, file = "ec.txt", sep = " ", row.names = FALSE, quote=FALSE)


test <- read_csv("Desktop/contigs/prodigal/reCOGnizer_results/14159.faa/test")

test <- test %>%
  dplyr::filter(date=="t2"|date=="t8") %>%
  group_by(`EC number`, date) %>%
  dplyr::mutate(tot=sum(norm_mapped_wa)) %>%
  dplyr::select(species,genus,date,tot,`EC number`) %>%
  distinct() %>%
  pivot_wider(names_from=date,values_from=tot) 

test$mg=1

unique(test$genus)

write_tsv(test,"Desktop/contigs/prodigal/reCOGnizer_results/14159.faa/test.tsv")
as.character(unique(test$`EC number`))