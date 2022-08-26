library(readr)
library(splitstackshape)
library(dplyr)



pig.id.basedir = "/Users/dgaio/Desktop"

for (pig.id in list.files(pig.id.basedir, pattern = "counts")[5]) {
  pig.id.dir = file.path(pig.id.basedir, pig.id)
  
  contig_counts <- read_table(pig.id.dir, 
                  col_names = FALSE)
  contig_counts$sample <- pig.id
  contig_counts$sample <- gsub("counts_",replacement = "",x=contig_counts$sample)
  
  contig_counts <- cSplit(contig_counts, splitCols = "sample","_")
  contig_counts <- cSplit(contig_counts, splitCols = "sample_2","-")
  
  contig_counts$sample_2_3 <- gsub(".bam",replacement = "",x=contig_counts$sample_2_3)
  
  

  #if pig/sow
  contig_counts <- contig_counts %>%
    select(sample_1,sample_2_2,sample_2_3,X1,X2,X3,X4)
  colnames(contig_counts) <- c ("pig", "date","duplicate", "contig","contigLen","mapped","unmapped")
  print(head(contig_counts))

  
  # # if pos ctrl
  # contig_counts <- contig_counts %>%
  #   select(sample_1,sample_2_2,sample_2_3,X1,X2,X3,X4)
  # contig_counts$sample_2_3 <- gsub("0",replacement = "",x=contig_counts$sample_2_3)
  # contig_counts$pig <- paste0(contig_counts$sample_1,"_R",contig_counts$sample_2_3)
  # contig_counts$duplicate <- "01" # neg controls don't have duplicates, but this col is added to match pig/sow samples
  # contig_counts <- contig_counts %>%
  #   select(pig,sample_2_2,duplicate,X1,X2,X3,X4)
  # colnames(contig_counts) <- c ("pig", "date", "duplicate","contig","contigLen","mapped","unmapped")
  # print(head(contig_counts))
  
  # # if neg ctrl
  # contig_counts$sample_2_2 <- gsub(".bam",replacement = "",x=contig_counts$sample_2_2)
  # 
  # contig_counts$pig <- paste0(contig_counts$sample_2_1,"_R",contig_counts$sample_2_2)
  # contig_counts$date <- "tNONE"
  # 
  # contig_counts <- contig_counts %>%
  #   select(pig,date,X1,X2,X3,X4)
  # contig_counts$duplicate <- "01" # neg controls don't have duplicates, but this col is added to match pig/sow samples
  # colnames(contig_counts) <- c ("pig", "date", "contig","contigLen","mapped","unmapped","duplicate")
  # contig_counts <- contig_counts %>%
  #   select(pig,date,duplicate,contig,contigLen,mapped,unmapped)
  # print(head(contig_counts))
  
  
}


# take a positive control as an example, or a sow sample (because only one time point)
# and look at the bias: longer contig ==> more mapped reads 
# and look at the removal of this bias with the norMalization 
test <- contig_counts %>%
  select(contig,mapped,contigLen)

test <- test %>%
  dplyr::filter(contigLen>1) %>% # this just removes the one row that samtools created to show unmapped reads (to any contig)
  mutate(new_bin = cut(contigLen, breaks=8))

test %>% group_by(new_bin) %>% 
  summarise(min_Len=min(contigLen),
            max_Len=max(contigLen),
            min_map=min(mapped),
            max_map=max(mapped))

# before normalization
test %>% 
  ggplot(., aes(x=new_bin,y=mapped))+
  geom_boxplot()

# after normalization
test %>%
  mutate(norm_mapped=mapped/contigLen) %>% 
  ggplot(., aes(x=new_bin,y=norm_mapped))+
  geom_boxplot()


# Proceed to normalize count data:
# 1. first add pseudocount (+1) 
# 2. second, normalize: norm_count=mapped/contigLen




library(stringr)
a <- contig_counts
a$date <- str_replace(a$date,"(\\d{2})(\\d{2})(\\d{2})$","\\1-\\2-\\3")
a$date <- paste0(20,a$date)
a$date <- as.Date(a$date)



reformat_dates_fun(a)



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
  my_df1$date <- plyr::mapvalues(as.character(my_df1$date), from, to)
  
  return(my_df1)
}











# merge all? 









