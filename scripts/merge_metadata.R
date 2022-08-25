
# This script takes in either a weighted contigs file 
# or a bins file (which comes from a weighted contigs file)
# and merges metadata to it. 

## STEP 1. Metadata files are read in. 
## STEP 2. Within a loop (iterates through all subjects), metadata is merged. 

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
source_dir = "/shared/homes/152324/metapigs_function/source_data/" # should contain: cohorts.xlsx, PigTrial_GrowthWtsGE.hlsx.xlsx, weights.csv, weights_final.csv

#test here: 
#pig.id.basedir = "/Users/dgaio/cloudstor/Gaio/out_new_test" 
#source_dir = "/Users/dgaio/cloudstor/Gaio/github/metapigs_function/source_data/" # should contain: cohorts.xlsx, PigTrial_GrowthWtsGE.hlsx.xlsx, weights.csv, weights_final.csv



########################################################
######### STEP 1 #######################################
########################################################

# Read in metadata. 
# read in cohorts file: 
cohorts <- read_excel(paste0(source_dir,"cohorts.xlsx"), 
                      col_names = c("cohort","pig"), skip = 1)

# read in pig details file & parse it:
pig_details <- read_excel(paste0(source_dir,"PigTrial_GrowthWtsGE.hlsx.xlsx"), 
                          sheet = "Piglet details", skip=1, 
                          col_types = c("text", "text", "skip", 
                                        "text", "text","date", "skip", 
                                        "text", "text", "text"),
                          col_names = c("pig","tattoo","crate","nurse_sow",
                                        "BIRTH_DAY","breed","sire","maternal_sow"))
pig_details$pig <- gsub("G","", pig_details$pig)
pig_details$pig <- gsub("T","", pig_details$pig)

# read in pig weight data & parse it:
weights <- read_csv(paste0(source_dir,"weights.csv"),
                    col_types = cols(Pig = col_character()))

weights_final <- read_csv(paste0(source_dir,"weights_final.csv"), 
                          col_types = cols(Pig = col_character()))

colnames(weights) <- c("room","pen","pig","t0","t2","t4","t6","t8")
colnames(weights_final) <- c("room","pen","pig","date","weight")

from = c("6-Mar","7-Mar","8-Mar","9-Mar","10-Mar")
to = c("t10","t10","t10","t10","t10")
weights_final$date <- plyr::mapvalues(as.character(weights_final$date), from, to) 
weights_final <- weights_final %>%
  dplyr::select(room,pen,pig,date,weight)

weights <- weights %>% dplyr::select(room,pen,pig,t0,t2,t4,t6,t8) %>%
  pivot_longer(., cols = c(t0,t2,t4,t6,t8), names_to = "date", values_to = "weight")

weights <- rbind(weights,weights_final) 


########################################################
######### STEP 2 #######################################
########################################################



# no loop, read in file and merge metadata 

# merge 1
f1 <- merge.data.frame(dff,cohorts) %>% 
  select("cohort", everything())

# merge 2
f2 <- left_join(f1,pig_details)

# merge 3
f3 <- left_join(f2,weights)

# reorder columns 
f4 <- f3 %>%
  select(cohort,pig,date,room,pen,weight,
         breed,BIRTH_DAY,crate,maternal_sow,nurse_sow,sire,
         bin,contigName,contigLen,value)


fwrite(
  x = f4,
  file = file.path(pig.id.dir, ##########),
  row.names=FALSE
)



