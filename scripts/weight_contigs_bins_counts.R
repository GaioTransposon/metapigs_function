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


pig.id.basedir = "/shared/homes/152324/contigs"
source_dir = "/shared/homes/152324/metapigs_function/source_data/" # should contain: Cdb.csv

# test on local 
#pig.id.basedir = "/Users/dgaio/cloudstor/Gaio/contigs"
#source_dir = "/Users/dgaio/cloudstor/Gaio/github/metapigs_function/source_data/" # should contain: Cdb.csv

#######
# Read in clustering info:

# clustering input file
Cdb = read.csv(paste0(source_dir,"Cdb.csv"))

# split genome column into pig id and bin
Cdb2 <- data.frame(Cdb, str_split_fixed(Cdb$genome, "_", 2))

# rename new columns
colnames(Cdb2)[colnames(Cdb2)=="X1"] <- "pig"
colnames(Cdb2)[colnames(Cdb2)=="X2"] <- "bin"

#subset df: pig, bin, secondary_cluster
Cdb3 <- Cdb2[c("pig", "bin", "secondary_cluster")]
#######

#######
# Function to get stats 
give_stats <- function(df) {
  
  type <- deparse(substitute(df))
  
  p <- unique(df$pig)
  d <- unique(df$date)
  
  df_out <- df %>%
    group_by(pig,date) %>%
    dplyr::summarise(pig=p, 
                     date=d, 
                     type=unique(type), 
                     mapped=sum(mapped),
                     unmapped=sum(unmapped),
                     contig_n=n(),
                     contig_len_tot=sum(contigLen),
                     contig_len_mean=mean(contigLen),
                     contig_len_med=median(contigLen))
  return(df_out)
}

# prep empty lists
my_list<-list() # for stats 
my_list2<-list() # for bins 

# open each file, get stats, weight mapped reads by contig length, save 
for (pig.id in list.files(pig.id.basedir, pattern = "parsed")) {
  pig.id.dir = file.path(pig.id.basedir, pig.id)

  counts_parsed <- read_csv(pig.id.dir,
                            col_types = cols(pig = col_character()))
  
  ##### stats
  # what proportion of contigs got binned? 
  no_contig_no_bin <- counts_parsed %>% 
    dplyr::filter(contig=="*")

  contig_only <- counts_parsed %>% 
    dplyr::filter(!contig=="*") %>%
    dplyr::filter(bin=="no_bin")
  
  contig_and_bin <- counts_parsed %>% 
    dplyr::filter(!contig=="*") %>%
    dplyr::filter(!bin=="no_bin")
  
  # if these sub-dataframes are not overlapping, then True: 
  NROW(no_contig_no_bin) + NROW(contig_only) + NROW(contig_and_bin) ==  NROW(counts_parsed)
  
  contig_stats <- rbind(as.data.frame(give_stats(contig_and_bin)),
        as.data.frame(give_stats(contig_only)),
        as.data.frame(give_stats(no_contig_no_bin))) %>%
    pivot_longer(cols=c(mapped,unmapped))
  
  my_list[[pig.id]] <- contig_stats
  #####
  
  ##### contigs 
  # normalize mapped read by contig length (POSSIBILITY: add pseudocopunt - discuss with Christian)
  weighted_contigs_counts <- counts_parsed %>% 
    dplyr::mutate(
      # mapped reads normalized by contig length
      mapped_wa=mapped/contigLen) 
  
  fwrite(weighted_contigs_counts, paste0(pig.id.dir, "_weighted_contigs.csv"))
  #####
  
  ##### bins 
  weighted_bins <- weighted_contigs_counts %>%
    dplyr::filter(!bin=="no_bin") %>%
    group_by(pig,date,bin) %>%
    dplyr::summarise(mapped_wa_bin=sum(mapped_wa))

  weighted_bins_clustered <- left_join(weighted_bins, Cdb3)
  
  my_list2[[pig.id]] <- weighted_bins_clustered
  #####
  
  
}


#######
# Stats save: 
contig_stats <- do.call("rbind", my_list)

contig_stats %>% 
  ggplot(., aes(x=type,y=value, fill=name), color=name) +
  geom_bar(stat="identity")

unique(contig_stats$date)

fwrite(contig_stats, paste0(dirname(pig.id.dir), "/contig_stats.csv"))
#######

#######
# Join weighted bins dataframes: 
all_bins <- do.call("rbind", my_list2)

fwrite(all_bins, paste0(dirname(pig.id.dir), "/all_wa_bins_counts.csv"))
#######


# # ######################################################
# # ######################################################
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
# weighted_contigs_counts %>%
#   dplyr::mutate(new_bin = cut(contigLen, breaks=8)) %>%
#   ggplot(., aes(x=new_bin,y=mapped_wa))+
#   geom_boxplot()
# 
# # ######################################################
# # ######################################################







