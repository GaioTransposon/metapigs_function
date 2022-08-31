library(readr)
library(dplyr)
library(ggplot2)
library(data.table)


#pig.id.basedir = "/shared/homes/152324/contigs/"

# test
pig.id.basedir = "/Users/dgaio/Desktop"


# open each file, add pseudocount, and normalize number of mapped reads by contig length 
for (pig.id in list.files(pig.id.basedir, pattern = "parsed")) {
  pig.id.dir = file.path(pig.id.basedir, pig.id)
  
  counts_parsed <- read_csv(pig.id.dir,
                            col_types = cols(pig = col_character()))
  
  # pseudocount addition
  counts_parsed$mapped <- counts_parsed$mapped+0.1 
  
  # weight contig counts 
  counts_weighted <- counts_parsed %>% 
    dplyr::filter(!contig=="*") %>% # remove row that contains unmapped reads data (reads that have not been mapped to any contig)
    dplyr::mutate(wa_mapped=mapped/contigLen) %>% # weighting mapped counts by contig length 
    dplyr::select(pig,date,contig,contigLen,wa_mapped)
  
  # save 
  fwrite(counts_weighted, paste0(pig.id.basedir, "/weighted_counts"))
  
}



# ######################################################
# ######################################################
# 
# # # Snippet of code to show why we need to weight the counts by the contig length:
# # # take a positive control as an example, or a sow sample (because only one time point
# # # so contigs with zero counts will not show up)
# # # and look at the bias: longer contig ==> more mapped reads
# # # and look at the removal of this bias with the normalization
# 
# # before normalization
# counts_parsed %>%
#   dplyr::mutate(new_bin = cut(contigLen, breaks=8)) %>%
#   ggplot(., aes(x=new_bin,y=mapped))+
#   geom_boxplot()
# 
# # after normalization
# counts_weighted %>%
#   dplyr::mutate(new_bin = cut(contigLen, breaks=8)) %>%
#   ggplot(., aes(x=new_bin,y=wa_mapped))+
#   geom_boxplot()
# 
# ######################################################
# ######################################################

