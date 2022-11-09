library(readr)
library(dplyr)
library(ggplot2)
library(plotly)
library(gapminder)

bin_depths <- read_csv("~/cloudstor/Gaio/github/metapigs_function/middle_dir/all_bins_depths.csv")
bin_counts <- read_csv("~/cloudstor/Gaio/github/metapigs_function/middle_dir/all_wa_bins_counts.csv")


View(bin_depths)
View(bin_counts)


anti <- anti_join(bin_depths,bin_counts)
View(anti)
unique(anti$pig) # bin counts include the positive control samples. bin depths don't have the pos control samples. 

bin_counts2 <- bin_counts[ ! bin_counts$pig %in% unique(anti$pig), ]
NROW(bin_counts2)

NROW(bin_depths)
NROW(bin_counts)
df <- merge(bin_counts,bin_depths, by = c("pig","date","bin","secondary_cluster")) %>% dplyr::filter(!date=="tNONE") #%>% dplyr::filter(pig=="29961")
NROW(df)
View(df)



p <- df %>%
  ggplot( aes(log(mapped_wa_bin), 
              log(value))) +
  geom_point(size=1) +
  theme_bw()

ggplotly(p)



