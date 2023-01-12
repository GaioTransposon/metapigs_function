library(readr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(EnvStats)
library(purrr)

df <- read_csv("/Users/dgaio/Desktop/contigs/prodigal/reCOGnizer_results/14159.faa/14159_reCOGnizer_results_eval_filtered_final.csv")

# reorder dates 
df$date  = factor(df$date, levels=c("t0",
                                    "t1", 
                                    "t2",
                                    "t3",
                                    "t4",
                                    "t5",
                                    "t6",
                                    "t7",
                                    "t8",
                                    "t9",
                                    "t10"))

df0 <- df %>%
  dplyr::filter(date=="t2"|date=="t8")

unique(df0$`Functional category`)

df1 <- df0 %>% 
  dplyr::filter(`Functional category`=="Carbohydrate transport and metabolism")
View(df1)
fam <- df1 %>%
  group_by(`Protein description`) %>%
  tally() %>%
  arrange(desc(n))

# for each unique protein (based on Protein description), check which ones are significantly diff (ab) between dates 
out <- split( df1 , f = df1$`Protein description` )
NROW(out)


res <- map(out, ~kruskal.test(norm_mapped_wa ~ date, data = .)) 
pval <- sapply(res, "[", "p.value")
class(pval)

# to df
y <- stack(pval)
y <- y %>%
  dplyr::filter(values<0.05)

y$ind=gsub(".p.value*","",y$ind)

y<-as.list(y$ind)
class(y)

NROW(y)


# run through list of signif results
df1_sig <- subset(df1, (`Protein description` %in% y))


to_plot <- split( df1_sig , f = df1_sig$`Protein description` )
NROW(to_plot)


#to_plot <- to_plot[1:3]


#Pdf
pdf('/Users/dgaio/Desktop/contigs/Example_carb.pdf')
#Loop
for (i in 1:39){
  
  z <- to_plot[30]

  z <- do.call(rbind.data.frame, z)

  p1 <- z %>%
    ggplot(., aes(x=date,y=log(norm_mapped_wa))) +
    geom_boxplot()+
    ggtitle(z$`Protein description`)+
    stat_n_text(vjust=-1)
  
  View(z)

  p2 <- z %>%
    ggplot(., aes(x=date,y=log(norm_mapped_wa))) +
    geom_boxplot()+
    facet_wrap(~bin)+
    stat_n_text(vjust=-1)

  both <- ggarrange(
    p1,p2,ncol=2)

  plot(p1)

}
dev.off()











z <- to_plot[47]

z <- do.call(rbind.data.frame, z)

p1 <- z %>%
  ggplot(., aes(x=date,y=log(norm_mapped_wa))) +
  geom_boxplot()+
  ggtitle(z$`Protein description`)+
  stat_n_text(vjust=-1) 

p2 <- z %>%
  ggplot(., aes(x=date,y=log(norm_mapped_wa))) +
  geom_boxplot()+
  ggtitle(z$species)+
  facet_wrap(~bin)+
  stat_n_text(vjust=-1) 

both <- ggarrange(
  p1,p2,ncol=2)

z$species
xx <- z %>%
  dplyr::select(species,norm_mapped_wa,bin)
View(xx)

################################################################################

# # Connection dbcan with new data - for discussion in manuscript: 
# # Take along all significant hits from dbcan (differentially represented CAZy between timepoints). 
# # Save enzyme IDs to a list. Search the recognizer_results files for items of this list 
# # and see what `Protein description` they correspond to. Do we see the same/similar trends? 
#   
# library(openxlsx)
# dbcan = read.xlsx("https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8549361/bin/mgen-7-0501-s002.xlsx",sheet=2)
# dbcan %>% dplyr::filter(p_value<=0.05) %>%
#   dplyr::select(enzID) %>%
#   distinct()
# # save as list
# 
# # go through reCOGnizer_results_eval_filtered.csv files and grep each item 
# # for example: 
# # cat reCOGnizer_results_eval_filtered.csv | grep -w "GH25" | cut -f 6 | head

