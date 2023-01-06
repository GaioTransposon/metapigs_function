library(readr)
library(dplyr)
library(ggplot2)
library(plotly)
library(gapminder)

bin_depths <- read_csv("~/github/metapigs_function/middle_dir/all_bins_depths.csv")
bin_counts <- read_csv("~/github/metapigs_function/middle_dir/all_wa_bins_counts.csv")
bin_counts_old <- read_csv("/Users/dgaio/github/metapigs_dry/middle_dir/no_reps_all.csv")


View(bin_depths)
View(bin_counts)
View(bin_counts_old)


anti <- anti_join(bin_counts_old,bin_counts)
View(anti)
unique(anti$pig) # bin counts include the positive control samples. bin depths don't have the pos control samples. 

bin_counts2 <- bin_counts[ ! bin_counts$pig %in% unique(anti$pig), ]
NROW(bin_counts2)

NROW(bin_depths)
NROW(bin_counts)
bin_counts_vs_bin_counts_old <- merge(bin_counts,bin_counts_old, by = c("pig","date","bin","secondary_cluster")) %>% dplyr::filter(!date=="tNONE") #%>% dplyr::filter(pig=="29961")
NROW(bin_counts_vs_bin_counts_old)
head(bin_counts_vs_bin_counts_old)

bin_counts_vs_bin_depths <- merge(bin_counts,bin_depths, by = c("pig","date","bin","secondary_cluster")) %>% dplyr::filter(!date=="tNONE") #%>% dplyr::filter(pig=="29961")
NROW(bin_counts_vs_bin_depths)
head(bin_counts_vs_bin_depths)

p <- bin_counts_vs_bin_counts_old %>%
  ggplot( aes(log(mapped_wa_bin), 
              log(value))) +
  geom_point(size=1) +
  theme_bw()+
  xlab("log(new_counts)")+
  ylab("log(old_counts)")+
  labs(title = "abundance calculation methods compared - read counts mapped to bins"
       #, caption = "Data source: ToothGrowth"
  )
p

p2 <- bin_counts_vs_bin_depths %>%
  ggplot( aes(log(mapped_wa_bin), 
              log(value))) +
  geom_point(size=1) +
  theme_bw()+
  xlab("log(new_counts)")+
  ylab("log(depths)")+
  labs(title = "abundance calculation methods compared - new counts vs depths"
       #, caption = "Data source: ToothGrowth"
  )
p2


cor.test(df$mapped_wa_bin, df$value, 
         method=c("pearson", "kendall", "spearman"))


