
# WHAT THIS SCRIPT DOES: 
# This script takes the weighted (wa) contigs data (wa_contigs.csv) 
# from each subject, parses it, and merges metadata to it. 

# BACKGROUND INFO:
# Depending on whether the "subject" is a
# pig/sow, or a positive or a negative control, it reformats the file 
# in a different way . This is because controls don't have sampling dates
# or they are not formatted like the pigs/sows. Another difference is that 
# positive and negative controls are technical replicates of each other, so 
# the formatting has to happen separately from pig/sow samples. 

# STEPS:
## STEP 1. three functions are written: 1 for pig/sow, 1 for POS, 1 for NEG controls. \
# in either of these three functions three main things happen: 
##### - reformatting of dataframe from wide to long
##### - dates are reformatted
##### - samples are de-duplicated
## STEP 2. Within a loop (iterates through all subjects), either of the three functions is run. 
## STEP 3. Metadata files are read in. 
## STEP 4. Within a loop (iterates through all subjects), metadata is merged. 



library(readr)
library(reshape)
library(tidyr)
library(splitstackshape)
library(readxl)
library(stringr)
library(dplyr)
library(data.table)


# input dir
#pig.id.basedir = "/shared/homes/152324/out_new" # make sure you have permission to write in the "subjectDIR"/metabat folder
#source_dir = "/shared/homes/152324/metapigs_function/source_dir/" # should contain: cohorts.xlsx, PigTrial_GrowthWtsGE.hlsx.xlsx, weights.csv, weights_final.csv


#test here: 
pig.id.basedir = "/Users/dgaio/cloudstor/Gaio/out_new_test" 
source_dir = "/Users/dgaio/cloudstor/Gaio/github/metapigs_function/source_data/" # should contain: cohorts.xlsx, PigTrial_GrowthWtsGE.hlsx.xlsx, weights.csv, weights_final.csv




########################################################
######### STEP 1 #######################################
########################################################


# Function 1 
fun_pig_sow <- function (my_file) {
  
  wa_contigs <- read.table(
    file = my_file, 
    header = TRUE, 
    sep = ",", 
    row.names = NULL
  )
  
  # depth columns into one single column
  wa_contigs_long <- wa_contigs %>%
    pivot_longer(., cols=ends_with("bam"))
  
  wa_contigs_long <- cSplit(wa_contigs_long, splitCols = "name",".")
  
  # remove "X" upstream
  wa_contigs_long$name_1 <- substring(wa_contigs_long$name_1, 2)
  # reformat to date
  wa_contigs_long$name_1 <- as.Date(wa_contigs_long$name_1, format = "%y%m%d")
  
  
  # group and reformat dates 
  from = c("2017-01-30", # mothers (sows) samples
           "2017-01-31","2017-02-01",
           "2017-02-03",
           "2017-02-06","2017-02-07","2017-02-08",
           "2017-02-10",
           "2017-02-14",
           "2017-02-16","2017-02-17",
           "2017-02-21",
           "2017-02-24", 
           "2017-02-28",
           "2017-03-03",
           "2017-03-06","2017-03-07","2017-03-08","2017-03-09","2017-03-10")
  
  to = c("tM", # mothers (sows) samples
         "t0","t0",
         "t1",
         "t2","t2","t2",
         "t3",
         "t4",
         "t5","t5",
         "t6",
         "t7", 
         "t8",
         "t9",
         "t10","t10","t10","t10","t10")
  
  # replace collection dates (date format) with groups of collection dates (character format)
  wa_contigs_long$date <- plyr::mapvalues(as.character(wa_contigs_long$name_1), from, to)
  
  
  #de-duplicate: { motive: some piglets could not be sampled on a given sampling date
  # (or low quality sample was taken, so sample was taken again the next day),
  # they were re-sampled the morning after (before treatment was administered)
  # for this reason we are going to make groups of sampling dates as follows } 
  no_reps_wa_contigs <- wa_contigs_long %>%
    group_by(pig, bin, contigName, contigLen, date) %>%
    dplyr::summarise(value=mean(value))
  
}


# Function 2 

fun_neg_ctrls <- function(my_file) {
  
  wa_contigs <- read.table(
    file = my_file, 
    header = TRUE, 
    sep = ",", 
    row.names = NULL
  )
  
  # depth columns into one single column
  wa_contigs_long <- wa_contigs %>%
    pivot_longer(., cols=ends_with("bam"))
  
  wa_contigs_long <- cSplit(wa_contigs_long, splitCols = "name",".")
  
  wa_contigs_long$name_1 <- gsub("X","R", wa_contigs_long$name_1)
  wa_contigs_long$pig <- paste0(wa_contigs_long$pig,"_",wa_contigs_long$name_1)
  wa_contigs_long$date <- "tNONE"
  no_reps_wa_contigs <- wa_contigs_long %>%
    group_by(pig, bin, contigName, contigLen, date) %>%
    dplyr::summarise(value=mean(value))
}


# Function 3 

fun_pos_ctrls <- function (my_file) {
  
  wa_contigs <- read.table(
    file = my_file, 
    header = TRUE, 
    sep = ",", 
    row.names = NULL
  )
  
  # depth columns into one single column
  wa_contigs_long <- wa_contigs %>%
    pivot_longer(., cols=ends_with("bam"))
  
  wa_contigs_long <- cSplit(wa_contigs_long, splitCols = "name",".")
  
  # remove "X" upstream
  wa_contigs_long$name_1 <- substring(wa_contigs_long$name_1, 2)
  # reformat to date
  wa_contigs_long$name_1 <- as.Date(wa_contigs_long$name_1, format = "%y%m%d")
  
  # group and reformat dates 
  from = c("2017-08-14","2018-01-24") # dates of pos ctrls
  to = c("tNONE","tNONE") # dates of pos ctrls
  
  # replace collection dates (date format) with groups of collection dates (character format)
  wa_contigs_long$date <- plyr::mapvalues(as.character(wa_contigs_long$name_1), from, to)
  
  # create replicate name 
  wa_contigs_long$pig <- paste0(wa_contigs_long$pig,"_",paste0("R",wa_contigs_long$name_2))
  
  no_reps_wa_contigs <- wa_contigs_long %>%
    group_by(pig, bin, contigName, contigLen, date) %>%
    dplyr::summarise(value=mean(value))
  
}


########################################################
######### STEP 2 #######################################
########################################################


# Run eitehr of the three functions above, iterating through all subjects
for (pig.id in list.files(pig.id.basedir)) {
  pig.id.dir = file.path(pig.id.basedir, pig.id)
  
  # if the sample contains digits it's either a pig or a sow sample
  if (str_detect(basename(pig.id.dir), "[0-9]")) {
    my_file1 <- paste(pig.id.dir,"metabat/wa_contigs.csv", sep="/")
    print(paste0(my_file1,"   it's pig"))
    
    if (file_test("-f", my_file1)) {
      print(paste0(my_file1,"   exists"))
      df1 <- fun_pig_sow(my_file1)
      fwrite(
        x = df1,
        file = file.path(pig.id.dir, "metabat/wa_contigs_clean.csv"),
        row.names=FALSE
      )
    } else { print(paste0(my_file1,"   does not exist"))}
    
  }
  

  # if the sample doesn't contain digits it's either a positive or negative control
  if (str_detect(basename(pig.id.dir), "[0-9]")==FALSE) {

    ######
    # separate paths for negative and pos controls
    
    # neg controls
    if (substr(basename(pig.id.dir),1,3)=="Neg") {
      my_file2 <- paste(pig.id.dir,"metabat/wa_contigs.csv", sep="/")
      print(paste0(my_file2,"   it's neg"))
      
      if (file_test("-f", my_file2)) {
        print(paste0(my_file2,"   exists"))
        df2 <- fun_neg_ctrls(my_file2)
        fwrite(
          x = df2,
          file = file.path(pig.id.dir, "metabat/wa_contigs_clean.csv"),
          row.names=FALSE
        )
      } else { print(paste0(my_file2,"   does not exist"))}
      
    }
    
    # pos controls
    else {
      my_file3 <- paste(pig.id.dir,"metabat/wa_contigs.csv", sep="/")
      print(paste0(my_file3,"   it's pos"))
      
      if (file_test("-f", my_file3)) {
        print(paste0(my_file3,"   exists"))
        df3 <- fun_pos_ctrls(my_file3)
        
        fwrite(
          x = df3,
          file = file.path(pig.id.dir, "metabat/wa_contigs_clean.csv"),
          row.names=FALSE
        )
      } else { print(paste0(my_file3,"   does not exist"))}
        
        
      
    }
    ######
  }
}


########################################################
######### STEP 3 #######################################
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
######### STEP 4 #######################################
########################################################

# Loop over (clean) wa contig files to merge metadata 
for (pig.id in list.files(pig.id.basedir)) {
  pig.id.dir = file.path(pig.id.basedir, pig.id)
  f <- paste(pig.id.dir,"metabat/wa_contigs_clean.csv", sep="/")
  
  # if the file exists, merge it to other dataframes:
  if (file_test("-f", f)) {
    print(paste0(f,"   exists"))
    
    dff <- read_csv(f, 
                    col_types = cols(pig = col_character()))
    
    
    # merge 1
    f_cohorts <- merge.data.frame(dff,cohorts) %>% 
      select("cohort", everything())
    
    # merge 2
    f_cohorts_deets <- left_join(f_cohorts,pig_details)
    
    # merge 3
    f_cohorts_deets_weights <- left_join(f_cohorts_deets,weights)
    
    # reorder columns 
    final_wa_contigs <- f_cohorts_deets_weights %>%
      select(cohort,pig,date,room,pen,weight,
             breed,BIRTH_DAY,crate,maternal_sow,nurse_sow,sire,
             bin,contigName,contigLen,value)

    
    fwrite(
      x = final_wa_contigs,
      file = file.path(pig.id.dir, "metabat/wa_contigs_depths.csv"),
      row.names=FALSE
    )
    
  } else { print(paste0(f,"   does not exist"))}
}
  
