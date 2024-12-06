library(readr)
library(dplyr)
library(ggplot2)

# list and read-in chunks 
my_files <- list.files(path="/Users/dgaio/Desktop", pattern="chunk")

dfList <- list()
for (f in my_files) {
  chunk <- read_csv(paste0("Desktop/",f), col_types = cols(pig = col_character()))
  
  # normalize value by sample
  chunk <- chunk %>% 
    group_by(date) %>%
    mutate(norm=value/sum(value))
  
  dfList[[f]] <- chunk
}

# concatenate them 
df <- bind_rows(dfList, .id = "column_label")
head(df)

#re-order factor levels
df$date <- factor(df$date, levels=c('t0', 't2', 't6', 't8', 't10'))

df %>%
  group_by(EC_number) %>%
  tally() %>%
  arrange(desc(n))

df %>%
  filter(date=="t0"|date=="t2"|date=="t8"|date=="t10") %>%
  filter(EC_number=="7.1.2.2") %>%
  ggplot(., aes(x=date,y=log(value), color=pig)) +
  geom_point()
unique(df$CDD_ID)


