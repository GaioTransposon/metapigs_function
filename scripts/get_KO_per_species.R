
if (getwd()!="/Users/dgaio") {
  .libPaths( c( .libPaths(), "/shared/homes/152324/R/library") )
  .libPaths( c( .libPaths(), "/shared/homes/152324/miniconda3/envs/recognizer_env/lib/R/library") )
} else {
  message("Running on local")
}



library(readr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(stringr)
library(tidyr)
library(data.table)
library(tidyverse)


################################################################################
################################################################################


# PART 1: 
# get number of KOs per species, for each pathway. 
print("# PART 1: get number of KOs per species, for each pathway.")


# paths to use: 
my_path <- "~/contigs/prodigal/eggnogg/KEGG/"
print(my_path)

these_files <- list.files(path = my_path, pattern="^all")
these_files <- these_files

df_list=list()
for (f in these_files) {
  
  f_path=paste0(my_path,f)
  df <- read_csv(f_path, col_types = cols(
    species = col_character(),
    pig = col_character()))
  
  # go ahead if pathway is not empty (KOs were found)
  if (NROW(df) > 0) {
    
    NROW(df)
    descr <- unique(df_merged$pathway_description)
    print(f)
    print(descr)
    
    essential <- df_merged %>%
      dplyr::select(KO,species) %>%
      distinct()
    
    # if species is NA, assign no_bin
    essential$species <- essential$species %>% replace_na('no_bin')
    
    save <- essential %>%
      group_by(species) %>%
      tally()
    
    df_list[[descr]] <- save
    
  }
  
}




all_dfs <- rbindlist(df_list, idcol = "paths") 

filename <- paste0(my_path,"KO_species/","KO_to_species.csv")
fwrite(all_dfs, file=filename)
head(all_dfs)
NROW(all_dfs)


################################################################################
################################################################################


# PART 2: 
# get number of KOs per species, for each pathway, using only contigs found within nearly complete MAGs. 
print("# PART 2: get number of KOs per species, for each pathway, using only contigs found within nearly complete MAGs.")

# open df containing only nearly complete MAGs
checkm_all_nearly <- read_delim("~/github/metapigs_dry/middle_dir/checkm_all_nearly",
                                "\t", escape_double = FALSE, trim_ws = TRUE)
checkm_all_nearly <- dplyr::filter(checkm_all_nearly, !grepl("Completeness",Completeness))
checkm_all_nearly$Completeness <- as.numeric(checkm_all_nearly$Completeness)
checkm_all_nearly$Contamination <- as.numeric(checkm_all_nearly$Contamination)
max(checkm_all_nearly$Contamination)
min(checkm_all_nearly$Completeness)
colnames(checkm_all_nearly) <- c("pig","bin","Marker lineage","Completeness","Contamination","Taxonomy (contained)")

checkm_all_nearly <- checkm_all_nearly %>% dplyr::select(pig,bin)


# paths to use: 
my_path <- "~/contigs/prodigal/eggnogg/KEGG/"
print(my_path)

these_files <- list.files(path = my_path, pattern="^all")
these_files <- these_files

df_list=list()
for (f in these_files) {
  
  f_path=paste0(my_path,f)
  df <- read_csv(f_path, col_types = cols(
    species = col_character(),
    pig = col_character()))
  
  # go ahead if pathway is not empty (KOs were found)
  if (NROW(df) > 0) {
    
    NROW(df)
    NROW(checkm_all_nearly)
    
    df_merged <- inner_join(df,checkm_all_nearly, by=c("pig","bin"))
    NROW(df_merged)
    
    # it can happen that no bins are left after the inner_join:
    if (NROW(df_merged) > 0) {
      
      descr <- unique(df_merged$pathway_description)
      print(f)
      print(descr)
      
      essential <- df_merged %>%
        dplyr::select(KO,species) %>%
        distinct()
      
      # if species is NA, assign no_bin
      essential$species <- essential$species %>% replace_na('no_bin')
      
      save <- essential %>%
        group_by(species) %>%
        tally()
      
      df_list[[descr]] <- save
      
    }
    
  }
  
}



all_dfs <- rbindlist(df_list, idcol = "paths") 

filename <- paste0(my_path,"KO_species/","KO_to_species_bestMAGs.csv")
fwrite(all_dfs, file=filename)
head(all_dfs)
NROW(all_dfs)


################################################################################
################################################################################
