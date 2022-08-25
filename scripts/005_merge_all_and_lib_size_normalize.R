

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
source_dir = "/shared/homes/152324/metapigs_function/source_data/" 

#test here: 
#pig.id.basedir = "/Users/dgaio/cloudstor/Gaio/out_new_test" 
#source_dir = "/Users/dgaio/cloudstor/Gaio/github/metapigs_function/source_data/" 




########################################################
######### STEP 2 #######################################
########################################################

# Loop over desired files and merge into a single DF 
for (pig.id in list.files(pig.id.basedir)) {
  pig.id.dir = file.path(pig.id.basedir, pig.id)
  f <- paste(pig.id.dir,#####desired file, sep="/")
  
  # if the file exists, merge it to other dataframes:
  if (file_test("-f", f)) {
    print(paste0(f,"   exists"))
    
    dff <- read_csv(f, 
                    col_types = cols(pig = col_character()))
    
    

    
    
    fwrite(
      x = final_wa_contigs,
      file = file.path(pig.id.dir, #####),
      row.names=FALSE
    )
    
  } else { print(paste0(f,"   does not exist"))}
}
