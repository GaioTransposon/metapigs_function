
# .libPaths( c( .libPaths(), "/shared/homes/152324/R/library") )
# .libPaths( c( .libPaths(), "/shared/homes/152324/miniconda3/envs/recognizer_env/lib/R/library") )


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


CDD_list=list()
for (f in these_files) {
  
  f_path=paste0(my_path,f)
  df <- read_csv(f_path)
  
  fname_new=str_replace(f, ".csv","_to_plot.csv")
  fnew_path=paste0(my_path,fname_new)
  
  to_plot <- split( df , f = df$`CDD ID` )
  
  for (i in seq(from = 1, to = NROW(to_plot), by = 1)) {
    
    CDD <- names(to_plot[i])
    pv <- to_plot[[i]] %>% dplyr::select(pvalue) %>% distinct()
    z <- to_plot[[i]] %>%
      group_by(date) %>%
      summarise_at(vars(norm_mapped_wa),
                   list(min=min, Q1=~quantile(., probs = 0.25),
                        median=median, Q3=~quantile(., probs = 0.75),
                        max=max))  
    t_before_med <- z[1,4]
    t_after_Q3 <- z[2,5]
    t_after_Q1 <- z[2,3]
    
    if (t_before_med >= t_after_Q3 || t_before_med <= t_after_Q1) {
      
      print("appending")
      CDD_list <- append(CDD_list,CDD)
      
    }
    
    else {
      
      print("not appended")

      
    }
    
  }
}
  
  

for (f in these_files) {
  
  f_path=paste0(my_path,f)
  df <- read_csv(f_path)

  # subset to proteins present in list made above: 
  df <- df[df$`CDD ID` %in% CDD_list,]
  
  fname_new=str_replace(f, ".csv","_major_shift.csv")
  fnew_path=paste0(my_path,fname_new)
  fwrite(df, file=fnew_path, sep = ",")
  
  if (NROW(df)!=0) {
    
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
      pdf_name=paste0(fnew_path,counter,'_subset_local.pdf')
      
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

  else {

    print("not printing a pdf")

  }

}


df %>%
  dplyr::filter(`CDD ID`=="CDD:227155") %>%
  group_by(date) %>%
  summarise_at(vars(norm_mapped_wa),
               list(min=min, Q1=~quantile(., probs = 0.25),
                    median=median, Q3=~quantile(., probs = 0.75),
                    max=max))  

library(agrmt)
a %>%
  dplyr::filter(`CDD ID`=="CDD:227155") %>%
  ggplot(., aes(x=date,y=log(norm_mapped_wa))) +
  geom_boxplot()+
  stat_n_text(vjust=-1)+
  theme(title = element_text(size=5))

a <- df %>%
  dplyr::filter(`CDD ID`=="CDD:227155") %>%
  group_by(date) %>%
  dplyr::mutate(norm_mapped_wa=norm_mapped_wa+minnz(norm_mapped_wa)/10)
View(a)






