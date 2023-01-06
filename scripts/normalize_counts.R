





# grep make sublists of files based on pig id 

# loop through sublists: 
#   for each file: 
#      normalize mapped_wa by date
#      concatenate into file, save based on pig id 



library(readr)
library(tidyverse)
library(dplyr)

# set directory: Desktop (or HPC) - contigs 
#pig.id.basedir = "~/contigs"

# local
pig.id.basedir = "~/Desktop/contigs"
these_files <- list.files(pig.id.basedir, pattern = "Counts_parsed_weighted_contigs")

names(these_files)


l1 <- list()
l2 <- list()

for (filename in these_files) {
  
  filenamesplit <- str_split(filename, "_",simplify = TRUE)
  pigID <- filenamesplit[1, 2]
  

  l1 <- append(l1,pigID)
  l2 <- append(l2,filename)
  
  
  my_path_to_file = file.path(pig.id.basedir, filename)
  df <- read_csv(my_path_to_file,
                 col_types = cols(pig = col_character()))
  
  my_sum <- sum(df$norm_mapped_wa)
  df_norm <- df %>% 
    dplyr::filter(!contig=="*") %>%
    group_by(pig,date,contig,bin) %>%
    dplyr::mutate(norm_mapped_wa=mapped_wa/sum(mapped_wa))
  
  head(df_norm)
  sum(df_norm$norm_mapped_wa)
  
  which(is.na(df$mapped_wa))
  View(df)
    

}


amap <- mapply(c,l1, l2, SIMPLIFY = FALSE)


for (i in 1:length(l1)){
  
  print (c(l1[i],l2[i]))
}







