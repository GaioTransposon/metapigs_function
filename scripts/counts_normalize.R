



# Snippet of code to show why we need to weight the counts by the contig length:

# # take a positive control as an example, or a sow sample (because only one time point)
# # and look at the bias: longer contig ==> more mapped reads 
# # and look at the removal of this bias with the normalization 
# test <- contig_counts %>%
#   select(contig,mapped,contigLen)
# 
# test <- test %>%
#   dplyr::filter(contigLen>1) %>% # this just removes the one row that samtools created to show unmapped reads (to any contig)
#   dplyr::mutate(new_bin = cut(contigLen, breaks=8))
# 
# test %>% group_by(new_bin) %>% 
#   summarise(min_Len=min(contigLen),
#             max_Len=max(contigLen),
#             min_map=min(mapped),
#             max_map=max(mapped))
# 
# # before normalization
# test %>% 
#   ggplot(., aes(x=new_bin,y=mapped))+
#   geom_boxplot()
# 
# # after normalization
# test %>%
#   mutate(norm_mapped=mapped/contigLen) %>% 
#   ggplot(., aes(x=new_bin,y=norm_mapped))+
#   geom_boxplot()












# open each and concatenate to existing dataframe (so each file won t stay in memory)
contig_counts %>% 
  dplyr::filter(!contig=="*") %>% # remove row that contains unmapped reads data (reads that have not been mapped to any contig)
  dplyr::select(pig,date,duplicate,contig,contigLen,mapped)


# open large file and run the following processes: 
# dedup (mean), pseudocount addition, and weighting contig counts 
contig_counts %>% 
  dplyr::filter(!contig=="*") %>%
  group_by(pig,date,contig) %>%
  dplyr::summarise(mapped=mean(mapped)+1, # deduplicating & adding pseudocount to all contigs 
                   contigLen=mean(contigLen), # deduplicating 
                   wa_mapped=mapped/contigLen) # weighting mapped counts by contig length 
fwrite(contig_counts, "weighted_counts")

# normalize based on lib size 
fwrite(weighted_counts_lib_size_normalized)




