library(readr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(EnvStats)
library(purrr)
library(data.table)

start_time <- Sys.time()

df <- read_csv("/Users/dgaio/Desktop/contigs/prodigal/reCOGnizer_results/14159.faa/14159_reCOGnizer_results_eval_filtered_final.csv")

# # reorder dates 
# df$date  = factor(df$date, levels=c("t0",
#                                     "t1", 
#                                     "t2",
#                                     "t3",
#                                     "t4",
#                                     "t5",
#                                     "t6",
#                                     "t7",
#                                     "t8",
#                                     "t9",
#                                     "t10"))

df0 <- df %>%
  dplyr::filter(date=="t2"|date=="t8")

df1 <- df0 #%>% 
  #dplyr::filter(`Functional category`=="Carbohydrate transport and metabolism")

# fam <- df1 %>%
#   group_by(`DB DB`) %>%
#   tally() %>%
#   arrange(desc(n))
# 
# test1=df1 %>%
#   dplyr::select(`Protein DB`)
# NROW(which(is.na(test1)))
# test2=df1 %>%
#   dplyr::select(`DB DB`)
# NROW(which(is.na(test2)))

# for each unique protein (based on Protein DB), check which ones are significantly diff (ab) between dates 
out <- split( df1 , f = df1$`DB description` )

#res <- map(out, ~kruskal.test(norm_mapped_wa ~ date, data = .)) 
res <- map(out, ~wilcox.test(norm_mapped_wa ~ date, data = .)) 

pval <- sapply(res, "[", "p.value")

# to df
y <- stack(pval)
y <- y %>%
  dplyr::filter(values<0.05)

y$ind=gsub(".p.value*","",y$ind)

# y<-as.list(y$ind)
# class(y)
# 
NROW(y)
# y
View(y)


# save results 
fwrite(y,"/Users/dgaio/Desktop/contigs/prodigal/reCOGnizer_results/14159.faa/14159_reCOGnizer_hits_R.tsv", sep = "\t")

end_time <- Sys.time()
end_time-start_time



z <- as.data.frame(X14159_reCOGnizer_hits_py)
colnames(z) <- c("ind","values")
z <- z %>% dplyr::select(values,ind)

View(y)
View(z)

class(z$values)

zz <- anti_join(y,z,by = c("ind"))

NROW(zz)
View(zz)







z[1,]==y[1,]

y %>% dplyr::filter(ind=="")
z %>% dplyr::filter(ind=="PTS sugar transporter subunit IIC.")







# run through list of signif results
df1_sig <- subset(df1, (`Protein DB` %in% y))

to_plot <- split( df1_sig , f = df1_sig$`Protein DB` )
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
    ggtitle(z$`Protein DB`)+
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


View(df1)

test <- df1 %>%
  dplyr::filter(`DB description`=='2 3 bisphosphoglycerate independent phosphoglycerate mutase iPGM. The 2,3-diphosphoglycerate- independent phosphoglycerate mutase (iPGM) catalyzes the interconversion of 3-phosphoglycerate (3PGA) and 2-phosphoglycerate (2PGA). They are the predominant PGM in plants and some other bacteria, including endospore forming Gram-positive bacteria and their close relatives. The two steps catalysis is a phosphatase reaction removing the phosphate from 2- or 3-phosphoglycerate, generating an enzyme-bound phosphoserine intermediate, followed by a phosphotransferase reaction as the phosphate is transferred from the enzyme back to the glycerate moiety. The iPGM exists as a dimer, each monomer binding 2 magnesium atoms, which are essential for enzymatic activity.')

View(test)




z <- to_plot[47]

z <- do.call(rbind.data.frame, z)

p1 <- z %>%
  ggplot(., aes(x=date,y=log(norm_mapped_wa))) +
  geom_boxplot()+
  ggtitle(z$`Protein DB`)+
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
# # and see what `Protein DB` they correspond to. Do we see the same/similar trends? 
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






t2_t8 <- read_delim("~/Desktop/contigs/prodigal/reCOGnizer_results/t2_t8", 
                    delim = "\t", escape_double = FALSE, trim_ws = TRUE)



x <- t2_t8 %>%
  group_by(protein,shift) %>%
  tally()
View(x)

colnames(df1)
df1 %>%
  dplyr::filter(`DB description`=="Uncharacterized conserved protein related to pyruvate formate-lyase activating enzyme  [Function unknown].") %>%
  dplyr::select(`species`) %>%
  group_by(species) %>%
  tally()


t2_t8 %>%
  dplyr::filter(protein=="Uncharacterized conserved protein related to pyruvate formate-lyase activating enzyme  [Function unknown].") %>%
  group_by(shift) %>%
  tally() 


df1 %>%
  dplyr::filter(`DB description`=="Uncharacterized conserved protein related to pyruvate formate-lyase activating enzyme  [Function unknown].") %>%
  ggplot(., aes(x=date,y=log(norm_mapped_wa))) +
  geom_boxplot()+
  facet_wrap(~bin)+
  stat_n_text(vjust=-1) 


t2_t8 %>%
  dplyr::select(protein) %>%
  distinct()

df1 %>%
  dplyr::select(`DB description`) %>%
  distinct()
  





