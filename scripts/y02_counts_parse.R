library(readr)
library(splitstackshape)
library(dplyr)
library(ggplot2)
library(stringr)
library(data.table)
library(tidyverse)


pig.id.basedir = "/shared/homes/152324/contigs" 

# test
# pig.id.basedir = "/Users/dgaio/Desktop/contigs"

# read in contigs to bins associations
all_bins_to_contigs <- read_csv(paste0(pig.id.basedir,"/all_bins_to_contigs.csv"))


# reformat dates function 
reformat_dates_fun <- function (my_df) {
  
  my_df1 <- my_df
  # reformat to date
  my_df1$date <- as.Date(my_df1$date, format = "%y%m%d")
  
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
           "2017-03-06","2017-03-07","2017-03-08","2017-03-09","2017-03-10",
           "2018-01-24","2017-08-14") # pos controls protexin (D-Scour), ColiGuard, and Mock Community
  
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
         "t10","t10","t10","t10","t10",
         "tNONE","tNONE") # pos controls protexin (D-Scour) and ColiGuard
  
  # replace collection dates (date format) with groups of collection dates (character format)
  my_df1$date <- plyr::mapvalues(as.character(my_df1$date), from, to)
  
  return(my_df1)
}



# list bam files 
my_list <- list.files(pig.id.basedir, pattern = "counts")
my_list <- list.files(pig.id.basedir, pattern = ".bam")
# exclude the output
my_list <-my_list[!str_detect(my_list,pattern="parsed")]
# exclude the dups 
duplicates <- read_table(paste0(pig.id.basedir,"/duplicates"), col_names = FALSE)
discard <- duplicates %>%
  group_by(X1) %>%
  dplyr::arrange(desc(X3)) %>%
  top_n(-1)
my_list <- my_list[!(my_list %in% discard$X2)]





# parse: 
for (pig.id in my_list) {
  pig.id.dir = file.path(pig.id.basedir, pig.id)
  
  contig_counts <- read_table(pig.id.dir,
                  col_names = FALSE)
  contig_counts$sample <- pig.id
  contig_counts$sample <- gsub("counts_",replacement = "",x=contig_counts$sample)

  contig_counts <- cSplit(contig_counts, splitCols = "sample","_")
  contig_counts <- cSplit(contig_counts, splitCols = "sample_2","-")

  contig_counts$sample_2_3 <- gsub(".bam",replacement = "",x=contig_counts$sample_2_3)
  
  
  # if pos ctrl:
  if (substr(basename(pig.id.dir),8,10)=="Moc"|substr(basename(pig.id.dir),8,10)=="Pro"|substr(basename(pig.id.dir),8,10)=="Col") {
    
    # if pos ctrl
    print("it's a pos ctrl")
    contig_counts <- contig_counts %>%
      dplyr::select(sample_1,sample_2_2,sample_2_3,X1,X2,X3,X4)
    contig_counts$sample_2_3 <- gsub("0",replacement = "",x=contig_counts$sample_2_3)
    contig_counts$pig <- paste0(contig_counts$sample_1,"_R",contig_counts$sample_2_3)
    contig_counts <- contig_counts %>%
      dplyr::select(pig,sample_2_2,X1,X2,X3,X4)
    colnames(contig_counts) <- c ("pig", "date","contig","contigLen","mapped","unmapped")
    # date reformatting
    contig_counts$date <- str_replace(contig_counts$date,"(\\d{2})(\\d{2})(\\d{2})$","\\1-\\2-\\3")
    contig_counts$date <- paste0(20,contig_counts$date) # 20 upstream of year (e.g. "2017")
    contig_counts$date <- as.Date(contig_counts$date)
    contig_counts <- reformat_dates_fun(contig_counts)
    
    
    # 1. merge bin data 
    # first split into two columns, otherwise it won't merge with the bin data (replicates not given there)
    contig_counts <- cSplit(contig_counts, splitCols = "pig","_")
    colnames(contig_counts)[colnames(contig_counts)=="pig_1"] <- "pig"
    contig_counts_bins <- left_join(contig_counts,all_bins_to_contigs, by=c("contig","pig"))
    
    # 2. if contig has not been binned to any bin, explicitly attribute "no_bin" string
    contig_counts_bins$bin <- contig_counts_bins$bin %>% replace_na("no_bin")
    
    # 3. merge back together 
    contig_counts_bins$pig <- paste0(contig_counts_bins$pig,"_",contig_counts_bins$pig_2)
    
    contig_counts_bins <- contig_counts_bins %>%
      dplyr::select(pig, date, bin, contig, contigLen, mapped, unmapped)

    fwrite(
      x = contig_counts_bins,
      file = paste0(sub('\\.bam$', '', pig.id.dir),"_contig_Counts_parsed"),
      row.names=FALSE
    )
    
  }
  
  # if neg ctrl:
  if (substr(basename(pig.id.dir),8,10)=="Neg") {

    print("it's a neg ctrl")
    contig_counts$sample_2_2 <- gsub(".bam",replacement = "",x=contig_counts$sample_2_2)
    
    # merge bin data
    colnames(contig_counts)[colnames(contig_counts)=="sample_1"] <- "pig"
    colnames(contig_counts)[colnames(contig_counts)=="X1"] <- "contig"
    contig_counts_bins <- left_join(contig_counts,all_bins_to_contigs, by=c("contig","pig"))
    
    # if contig has not been binned to any bin, explicitly attribute "no_bin" string
    contig_counts_bins$bin <- contig_counts_bins$bin %>% replace_na("no_bin")
    contig_counts_bins$pig <- NULL
    contig_counts_bins$pig <- paste0(contig_counts_bins$sample_2_1,"_R",contig_counts_bins$sample_2_2)
    contig_counts_bins$date <- "tNONE"
    contig_counts_bins <- contig_counts_bins %>%
      dplyr::select(pig,date,contig,X2,X3,X4,bin)
    colnames(contig_counts_bins) <- c ("pig", "date", "contig","contigLen","mapped","unmapped","bin")
    print(head(contig_counts_bins))
    
    fwrite(
      x = contig_counts_bins,
      file = paste0(sub('\\.bam$', '', pig.id.dir),"_contig_Counts_parsed"),
      row.names=FALSE
    )

  }
  
  # if the substring contains some numbers it must be pig/sow:
  if (str_detect(substr(basename(pig.id.dir),8,14), "[0-9]")==TRUE) {
    
    print("it's a pig/sow")
    #if pig/sow
    contig_counts <- contig_counts %>%
      dplyr::select(sample_1,sample_2_2,X1,X2,X3,X4)
    colnames(contig_counts) <- c ("pig", "date", "contig","contigLen","mapped","unmapped")
    # date reformatting
    contig_counts$date <- str_replace(contig_counts$date,"(\\d{2})(\\d{2})(\\d{2})$","\\1-\\2-\\3")
    contig_counts$date <- paste0(20,contig_counts$date) # 20 upstream of year (e.g. "2017")
    contig_counts$date <- as.Date(contig_counts$date)
    contig_counts <- reformat_dates_fun(contig_counts)

    # merge bin data 
    contig_counts$pig <- as.character(contig_counts$pig)
    contig_counts_bins <- left_join(contig_counts,all_bins_to_contigs, by=c("contig","pig"))
    
    # if contig has not been binned to any bin, explicitly attribute "no_bin" string
    contig_counts_bins$bin <- contig_counts_bins$bin %>% replace_na("no_bin")
    
    print(head(contig_counts_bins))
    fwrite(
      x = contig_counts_bins,
      file = paste0(sub('\\.bam$', '', pig.id.dir),"_contig_Counts_parsed"),
      row.names=FALSE
    )
    
  }
  
}






