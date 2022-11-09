library(readr)
library(dplyr)
library(ggplot2)
library(readr)
library(splitstackshape)
library(dplyr)
library(ggplot2)
library(stringr)
library(data.table)
library(tidyverse)


#pig.id.basedir = "/shared/homes/152324/contigs/"

# test
pig.id.basedir = "/Users/dgaio/Desktop"





# open each file, add pseudocount, and normalize number of mapped reads by contig length 
for (pig.id in list.files(pig.id.basedir, pattern = "parsed")) {
  pig.id.dir = file.path(pig.id.basedir, pig.id)
  
  counts_parsed <- read_csv(pig.id.dir,
                            col_types = cols(pig = col_character()))
  
  
  ####
  # get some stats 
  
  # save unmapped info (these are counts not mapped to any contig)
  unmapped_to_any_contig <- counts_parsed %>% 
    dplyr::filter(contig=="*")
  
  # info contigs per bin 
  stats_contig_per_bin <- counts_parsed %>% 
    dplyr::filter(!contig=="*") %>% # remove "unmapped" to any contig 
    group_by(pig,date,bin) %>%
    summarise(no_contigs=n(),
      # number of contigs per bin 
      mean_contigLen=mean(contigLen),
      # contig length info
      median_contigLen=median(contigLen),
      sd_contigLen=sd(contigLen),
      min_contigLen=min(contigLen),
      max_contigLen=max(contigLen),
      # mapping info 
      mean_mapped=mean(mapped),
      mean_unmapped=mean(unmapped),
      median_mapped=median(mapped),
      median_unmapped=median(unmapped),
      sd_mapped=sd(mapped),
      sd_unmapped=sd(unmapped),
      min_mapped=min(mapped),
      max_mapped=max(mapped),
      min_unmapped=min(unmapped),
      max_unmapped=max(unmapped))
  
  # stats of all bins 
  stats_overall <- counts_parsed %>% 
    group_by(pig,date) %>%
    dplyr::filter(!contig=="*") %>% # remove "unmapped" to any contig 
    summarise(no_contigs=n(),
              # number of contigs per bin 
              mean_contigLen=mean(contigLen),
              # contig length info
              median_contigLen=median(contigLen),
              sd_contigLen=sd(contigLen),
              min_contigLen=min(contigLen),
              max_contigLen=max(contigLen),
              # mapping info 
              mean_mapped=mean(mapped),
              mean_unmapped=mean(unmapped),
              median_mapped=median(mapped),
              median_unmapped=median(unmapped),
              sd_mapped=sd(mapped),
              sd_unmapped=sd(unmapped),
              min_mapped=min(mapped),
              max_mapped=max(mapped),
              min_unmapped=min(unmapped),
              max_unmapped=max(unmapped))

  ####
  
  # these are all unique contigs 
  NROW(unique(counts_parsed$contig))==NROW(counts_parsed)
  
  counts_weighted <- counts_parsed %>% 
    dplyr::filter(!contig=="*") %>% # remove "unmapped" to any contig --> careful cause here you ignore unmapped reads (to any contig/bin)
    # we only need to weight the mapped counts by the contig length
    dplyr::mutate(wa_mapped=mapped/contigLen) 
  
  ###
  # for bins: 
  weighted_bins_counts <- counts_weighted %>% 
    group_by(pig,date,bin) %>%
    # and sum these by bin. No need for pseudocount addition because no bins end up having zero counts
    dplyr::summarise(wa_mapped_bin=sum(wa_mapped))
  # fuck! >85% is in no_bin!!!!
  
  # # visualize 
  # weighted_bins_counts %>%
  #   filter(wa_mapped_bin<6000) %>%
  #   ggplot(., aes(x=reorder(bin,wa_mapped_bin),y=wa_mapped_bin))+
  #   geom_point()
  ###
  
  ###
  # for contigs: 
  weighted_contigs_counts <- counts_weighted %>% 
    dplyr::mutate(
      # and add a pseudocount, whcih is the minimum non-zero mapped value divided by 1000
      wa_mapped_pseudo=wa_mapped+(min(wa_mapped[wa_mapped > 0])/1000))
  ###

  # save 
  # remove last n chars from file name 
  fwrite(weighted_contigs_counts, paste0(pig.id.dir, "_weighted_contig_counts"))
  fwrite(weighted_bins_counts, paste0(pig.id.dir, "_weighted_bin_counts"))
  
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

