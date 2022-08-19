
# WHAT THIS SCRIPT DOES: 
# This script takes the weighted (wa) contigs data with the metadata merged to it 
# (wa_contigs_depths.csv)
# from each subject, and summarises it into bin depths. 

# STEPS:
## STEP 1. contigs depths to bins depths 
## STEP 2. Reads in the clustering info file (output of dRep) & 
# merges clustering info to bin depths. 


library(readr)
library(reshape)
library(tidyr)
library(splitstackshape)
library(readxl)
library(stringr)
library(dplyr)
library(data.table)


# input dir
pig.id.basedir = "/shared/homes/152324/out_new" # make sure you have permission to write in the "subjectDIR"/metabat folder
source_dir = "/shared/homes/152324/metapigs_function/source_data/" # should contain: Cdb.csv

#test here: 
#pig.id.basedir = "/Users/dgaio/cloudstor/Gaio/out_new_test" 
#source_dir = "/Users/dgaio/cloudstor/Gaio/github/metapigs_function/source_data/" # should contain: Cdb.csv



########################################################
######### STEP 1 #######################################
########################################################


# WA contigs depths to bins depths 
for (pig.id in list.files(pig.id.basedir)[4:5]) {
  pig.id.dir = file.path(pig.id.basedir, pig.id, "metabat")
  
  f <- paste(pig.id.dir, "wa_contigs_clean.csv", sep="/")

  # if the file exists, open it
  if (file_test("-f", f)) {
    print(paste0(f,"   file exists"))
    
    df <- read.table(
      file = f, 
      header = TRUE, 
      sep = ",", 
      row.names = NULL
    )
    
    # # take the mean or median?
    # check <- df %>%
    #   group_by(pig,bin) %>%
    #   summarise(mea=mean(value),
    #             med=median(value))
    # summary(check$mea)
    # summary(check$med)
    # #mean better (otherwise you create a bunch of zeros)
    

    wa_bins <- df %>%
      group_by(pig,bin,date) %>% 
      summarise(value=mean(value),
                binLen=sum(contigLen),
                binLen_sd=sd(contigLen))
    
    # # little sanity check
    # check <- df %>%
    #   dplyr::filter(bin=="bins.1.fa") %>%
    #   dplyr::filter(date=="t0")
    # out <- wa_bins %>%
    #   dplyr::filter(bin=="bins.1.fa") %>%
    #   dplyr::filter(date=="t0")
    # out$value==mean(check$value)
    # out$binLen==sum(check$contigLen)
    
    fwrite(
      x = wa_bins,
      file = file.path(pig.id.dir, "wa_bins.csv"),
      row.names=FALSE)
    
  }
  
  else { print(paste0(f,"   file does not exist"))}
  
}

    
    





########################################################
######### STEP 2 #######################################
########################################################


# 1. clustering input file
Cdb = read.csv(paste0(source_dir,"Cdb.csv"))

# split genome column into pig id and bin
Cdb2 <- data.frame(Cdb, str_split_fixed(Cdb$genome, "_", 2))

# rename new columns
colnames(Cdb2)[colnames(Cdb2)=="X1"] <- "pig"
colnames(Cdb2)[colnames(Cdb2)=="X2"] <- "bin"

#subset df: pig, bin, secondary_cluster
Cdb3 <- Cdb2[c("pig", "bin", "secondary_cluster")]


# 6. merging bins clustering info to depths, creating a file in each
# pig dir containing the weighted depths per bin and its clustering info
for (pig.id in list.files(pig.id.basedir)) {
  pig.id.dir = file.path(pig.id.basedir, pig.id, "metabat")
  
  f <- paste(pig.id.dir, "wa_bins.csv", sep="/")
  
  # if the file exists, open it
  if (file_test("-f", f)) {
    print(paste0(f,"   file exists"))
    
    # read in wa_bins.csv
    wa_bins_df <- read.table(
      file = f, 
      header = TRUE, 
      sep = ",", 
      row.names = NULL
    )
    wa_bins_df$pig <- as.character(wa_bins_df$pig)

    # merge
    clu_wa_bins <- dplyr::left_join(wa_bins_df, Cdb3, by=c("pig","bin"))
    
    #move secondary_cluster column to first position of dataframe
    clustered_wa_bins <- clu_wa_bins %>% 
      dplyr::select(secondary_cluster, everything())
    
    fwrite(
      x = clustered_wa_bins,
      file = file.path(pig.id.dir, "clustered_wa_bins.csv"),
      row.names=FALSE)
    
    
  }
  
  else { print(paste0(f,"   file does not exist"))}
  
}



