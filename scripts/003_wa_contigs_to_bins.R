


# 1. clustering input file
Cdb = read.csv(
  file="/shared/homes/12705859/Cdb.csv"
)
#test here: 
# Cdb = read.csv(
#   file="/Users/dgaio/cloudstor/Gaio/github/metapigs_dry/middle_dir/Cdb.csv"
# )



#5. Length-normalized bins' depths (summed contigs' depths) 
for (pig.id in list.files(pig.id.basedir)) {
  pig.id.dir = file.path(pig.id.basedir, pig.id, "metabat")
  
  wa_contigs <- paste(pig.id.dir, "wa_contigs.csv", sep="/")
  
  df <- read.table(
    file = wa_contigs, 
    header = TRUE, 
    sep = ",", 
    row.names = NULL
  )
  
  # # take the mean or median? 
  # df %>%
  #   group_by(pig,bin) %>%
  #   summarise(mea=mean(X170131.01.bam),
  #             med=median(X170131.01.bam))
  # #mean better (otherwise you create a bunch of zeros)
  
  cols_with_bam <- df %>% 
    select(contains('.bam'))
  cols_with_bam <- colnames(cols_with_bam)
  
  wa_bins <- df %>%
    group_by(pig,bin) %>%
    summarise(across(cols_with_bam, mean),
              binLen=sum(contigLen))
  
  # # little sanity check
  # check <- df %>%
  #   filter(bin=="bins.1.fa") %>%
  #   select(X170131.01.bam,contigLen)
  # out <- wa_bins %>%
  #   filter(bin=="bins.1.fa") %>%
  #   select(X170131.01.bam,binLen)
  # out[1,2]==mean(check$X170131.01.bam)
  # sum(check$contigLen)==out$binLen
  
  fwrite(
    x = wa_bins,
    file = file.path(pig.id.dir, "wa_bins.csv"),
    row.names=FALSE)
  
}


# split genome column into pig id and bin
Cdb2 <- data.frame(Cdb, str_split_fixed(Cdb$genome, "_", 2))

# rename new columns
colnames(Cdb2)[colnames(Cdb2)=="X1"] <- "pig"
colnames(Cdb2)[colnames(Cdb2)=="X2"] <- "bin"

#subset df: pig, bin, secondary_cluster
Cdb3 <- Cdb2[c("pig", "bin", "secondary_cluster")]

#change to character pig and bin items in dataframes
Cdb4 <- Cdb3 %>%  
  dplyr::mutate(pig = as.character(pig))
Cdb5 <- Cdb4 %>%  
  dplyr::mutate(bin = as.character(bin))

# 6. merging bins clustering info to depths, creating a file in each
# pig dir containing the weighted depths per bin and its clustering info
for (pig.id in list.files(pig.id.basedir)) {
  pig.id.dir = file.path(pig.id.basedir, pig.id, "metabat")
  
  # read in wa_bins.csv
  wa_bins <- paste(pig.id.dir, "wa_bins.csv", sep="/")
  wa_bins_df <- read.table(
    file = wa_bins, 
    header = TRUE, 
    sep = ",", 
    row.names = NULL
  )
  wa_bins_df2 <- wa_bins_df %>%  
    dplyr::mutate(pig = as.character(pig))
  wa_bins_df3 <- wa_bins_df2 %>%  
    dplyr::mutate(bin = as.character(bin))
  
  # merge with Cdb3
  clu_wa_bins <- dplyr::left_join(wa_bins_df3, Cdb5, by=c("pig","bin"))
  
  #move secondary_cluster column to first position of dataframe
  clustered_wa_bins <- clu_wa_bins %>% 
    dplyr::select(secondary_cluster, everything())
  
  fwrite(
    x = clustered_wa_bins,
    file = file.path(pig.id.dir, "clustered_wa_bins.csv"),
    row.names=FALSE
  )
}




# fwrite file

# depths --> counts & write file 





