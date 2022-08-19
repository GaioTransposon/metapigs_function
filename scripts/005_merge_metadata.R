

## STEP 3. Metadata files are read in. 
## STEP 4. Within a loop (iterates through all subjects), metadata is merged. 




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
