library(readr)
rec <- read_csv("Desktop/contigs/prodigal/reCOGnizer_results/reCOGnizer_POI_METABOLISM_t2_vs_t8_significant_major_shift.csv")

a <- rec %>%
  group_by(pig, date,`CDD ID`) %>%
  dplyr::summarise(x=sum(norm_mapped_wa)) 

a$sample=paste0(a$date,"_",a$pig)
a$pig <- NULL
a$date <- NULL

b <- a %>%
  pivot_wider(names_from=`CDD ID`, values_from=x, values_fill = 0)

b <- as.data.frame(b)
rownames(b) <- b$sample
b$sample <-NULL
pca <- prcomp(b)
head(pca)

dates <- substr(rownames(b), start = 1, stop = 2)  
dates  = factor(dates, levels=c("t2","t8"))

# plot! 
p <- fviz_pca_ind(pca, 
                     geom.ind="point",
                     fill.ind = dates, 
                     pointshape = 21, 
                     pointsize = 2,
                  axes = c(1,2),  # PCAs to plot
                     habillage = dates,
                     #geom.ind = "point", # show points only (nbut not "text") 
                     col.ind = dates, # color by groups
                     palette = c("#00AFBB", "#FC4E07"),   #"#E7B800"
                     addEllipses = FALSE, # Concentration ellipses
                     title="")+
  # scale_color_manual(name="time point", 
  #                    values=rainbow(n = 2))+
  theme(legend.position="none")+
  guides(color = guide_legend(nrow = 1))




fviz_pca_biplot(pca,
                geom.ind="point",
                pointshape = 21, pointsize = 2,
                habillage = dates,
                alpha.var ="contrib",
                repel = TRUE,labelsize=4) +
  scale_color_manual(name="time point",
                     values=rainbow(n = 11))+
  ggtitle("")+
  theme(legend.position="right",
        panel.border = element_rect(colour = "black", fill=NA, size=1))




# together
together <- fviz_pca_biplot(p,
                            geom.ind="point",
                            pointshape = 21, pointsize = 2,
                            habillage = dates,
                            alpha.var ="contrib",
                            repel = TRUE,labelsize=4) +
  scale_color_manual(name="time point",
                     values=rainbow(n = 2))+
  ggtitle("")+
  theme(legend.position="right",
        panel.border = element_rect(colour = "black", fill=NA, size=1))

rec %>%
  dplyr::select(`Protein description`) %>%
  distinct()
  


