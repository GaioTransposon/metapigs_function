
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(patchwork)
# downloaded ComplexHeatmamp locally from: https://cran.r-project.org/src/contrib/Archive/rjson/
library(dendextend)
library(ComplexHeatmap)
library("cluster")
library(scales)
library(stringr)

dfs <- read.csv(file = '/Users/dgaio/contigs/prodigal/eggnogg/KEGG/heatmap.csv')
dfs <- dfs %>% dplyr::select(pathway_description,pathway,KO,t0,t2,t4,t6,t8,t10)
unique(dfs$pathway_description)

dfs <- subset(dfs, nchar(as.character(pathway_description)) <= 100)



give_heatmap <- function(mykeywords_selection){
  
  mykeywords_selection <- paste(mykeywords_selection, collapse='|')
  
  dfs_sub <- dfs[grepl(mykeywords_selection, dfs$pathway_description, ignore.case = TRUE), ]
  
  rownames(dfs_sub) <- paste0(dfs_sub$KO,'_',dfs_sub$pathway)
  mylabels <- dfs_sub$pathway_description
  
  dfs_sub$pathway_description <- NULL
  dfs_sub$pathway <- NULL
  dfs_sub$KO <- NULL
  
  m <- as.matrix(dfs_sub)
  
  col.order <- c("t0","t2","t4","t6","t8","t10") # "t0","t2","t4","t6","t8","t10"
  m <- m[ , col.order]
  
  m <- t(apply(m, 1, function(x) rescale(x, to=c(-1,1))))
  
  
  # split by a vector specifying rowgroups
  Heatmap(m, name = "avg_ab",
          split = mylabels, 
          row_names_gp = gpar(fontsize = 4), 
          cluster_row_slices = TRUE, 
          clustering_distance_rows = "euclidean",
          show_row_dend = FALSE,
          cluster_columns = FALSE, 
          width = unit(6, "cm"), 
          row_title_rot = 0, 
          column_names_rot = 0, gap = unit(0.05, "cm"),
          border = "black",
          row_title_gp = gpar(fontsize = 7), 
          cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(sprintf("%.1f", m[i, j]), x, y, gp = gpar(fontsize = 5))
          })
  
}



#abx <- c("antimicrobial", "antibiotic", "vancomycin", "CAMP", "drug resistance", "Lactam", "Antifolate")
#give_heatmap(abx)

# https://www.genome.jp/kegg-bin/show_organism?menu_type=pathway_maps&category=Bacteria

carb <- c("Glycolysis / Gluconeogenesis","Citrate cycle (TCA cycle)",
          "Pentose phosphate pathway","Pentose and glucuronate interconversions",
          "Fructose and mannose metabolism","Galactose metabolism",
          "Ascorbate and aldarate metabolism","Starch and sucrose metabolism",
          "Amino sugar and nucleotide sugar metabolism","Pyruvate metabolism",
          "Glyoxylate and dicarboxylate metabolism","Propanoate metabolism",
          "Butanoate metabolism","C-Branched dibasic acid metabolism",
          "Inositol phosphate metabolism")

ene <- c("Oxidative phosphorylation","Photosynthesis",
         "Photosynthesis - antenna proteins",
         "Carbon fixation in photosynthetic organisms",
         "Carbon fixation pathways in prokaryotes","Methane metabolism",
         "Nitrogen metabolism","Sulfur metabolism")

lipid <- c("Fatty acid biosynthesis","Fatty acid elongation",
           "Fatty acid degradation","Cutin, suberine and wax biosynthesis",
           "Steroid biosynthesis",
           "Primary bile acid biosynthesis","Secondary bile acid biosynthesis",
           "Steroid hormone biosynthesis","Glycerolipid metabolism",
           "Glycerophospholipid metabolism","Ether lipid metabolism",
           "Sphingolipid metabolism","Arachidonic acid metabolism",
           "Linoleic acid metabolism","alpha-Linolenic acid metabolism",
           "Biosynthesis of unsaturated fatty acids")


bile <- c("Primary bile acid biosynthesis","Secondary bile acid biosynthesis")

nuc <- c("Purine metabolism","Pyrimidine metabolism")

aa <- c("Alanine, aspartate and glutamate metabolism",
        "Glycine, serine and threonine metabolism","Cysteine and methionine metabolism",
        "Valine, leucine and isoleucine degradation",
        "Valine, leucine and isoleucine biosynthesis","Lysine biosynthesis",
        "Lysine degradation","Arginine biosynthesis",
        "Arginine and proline metabolism","Histidine metabolism","Tyrosine metabolism",
        "Phenylalanine metabolism","Tryptophan metabolism",
        "Phenylalanine, tyrosine and tryptophan biosynthesis")

oth_aa <- c("beta-Alanine metabolism","Taurine and hypotaurine metabolism",
            "Phosphonate and phosphinate metabolism", "Selenocompound metabolism",
            "Cyanoamino acid metabolism","D-Amino acid metabolism",
            "Glutathione metabolism")

glycan <- c("N-Glycan biosynthesis","Various types of N-glycan biosynthesis",
         "Mucin type O-glycan biosynthesis","Mannose type O-glycan biosynthesis",
         "Other types of O-glycan biosynthesis",
         "Glycosaminoglycan biosynthesis - chondroitin sulfate / dermatan sulfate",
         "Glycosaminoglycan biosynthesis - heparan sulfate / heparin",
         "Glycosaminoglycan biosynthesis - keratan sulfate",
         "Glycosaminoglycan degradation","Glycosylphosphatidylinositol (GPI)-anchor biosynthesis",
         "Glycosphingolipid biosynthesis - lacto and neolacto series",
         "Glycosphingolipid biosynthesis - globo and isoglobo series",
         "Glycosphingolipid biosynthesis - ganglio series",
         "Other glycan degradation","Lipopolysaccharide biosynthesis",
         "O-Antigen repeat unit biosynthesis","O-Antigen nucleotide sugar biosynthesis",
         "Peptidoglycan biosynthesis","Teichoic acid biosynthesis",
         "Lipoarabinomannan (LAM) biosynthesis",
         "Arabinogalactan biosynthesis - Mycobacterium",
         "Exopolysaccharide biosynthesis") 

co_vit <- c("Thiamine metabolism","Riboflavin metabolism","Vitamin B metabolism",
            "Nicotinate and nicotinamide metabolism","Pantothenate and CoA biosynthesis",
            "Biotin metabolism","Lipoic acid metabolism","Folate biosynthesis",
            "One carbon pool by folate","Retinol metabolism","Porphyrin metabolism",
            "Ubiquinone and other terpenoid-quinone biosynthesis")


terp <- c("Terpenoid backbone biosynthesis","Monoterpenoid biosynthesis",
          "Sesquiterpenoid and triterpenoid biosynthesis","Diterpenoid biosynthesis",
          "Carotenoid biosynthesis","Brassinosteroid biosynthesis",
          "Insect hormone biosynthesis","Zeatin biosynthesis","Limonene degradation",
          "Pinene, camphor and geraniol degradation","Type I polyketide structures",
          "Biosynthesis of 12-, 14- and 16-membered macrolides","Biosynthesis of ansamycins",
          "Biosynthesis of enediyne antibiotics","Biosynthesis of type II polyketide backbone",
          "Biosynthesis of type II polyketide products","Tetracycline biosynthesis",
          "Polyketide sugar unit biosynthesis","Nonribosomal peptide structures",
          "Biosynthesis of siderophore group nonribosomal peptides",
          "Biosynthesis of vancomycin group antibiotics")

sec_met <- c("Phenylpropanoid biosynthesis",
             "Stilbenoid, diarylheptanoid and gingerol biosynthesis","Flavonoid biosynthesis",
             "Flavone and flavonol biosynthesis","Anthocyanin biosynthesis",
             "Isoflavonoid biosynthesis","Degradation of flavonoids",
             "Indole alkaloid biosynthesis","Indole diterpene alkaloid biosynthesis",
             "Isoquinoline alkaloid biosynthesis",
             "Tropane, piperidine and pyridine alkaloid biosynthesis",
             "Biosynthesis of various alkaloids","Caffeine metabolism","Betalain biosynthesis",
             "Glucosinolate biosynthesis","Benzoxazinoid biosynthesis",
             "Penicillin and cephalosporin biosynthesis","Carbapenem biosynthesis",
             "Monobactam biosynthesis","Clavulanic acid biosynthesis",
             "Streptomycin biosynthesis","Neomycin, kanamycin and gentamicin biosynthesis",
             "Acarbose and validamycin biosynthesis","Novobiocin biosynthesis",
             "Staurosporine biosynthesis","Phenazine biosynthesis","Prodigiosin biosynthesis",
             "Aflatoxin biosynthesis","Biosynthesis of various antibiotics",
             "Biosynthesis of various plant secondary metabolites",
             "Biosynthesis of various other secondary metabolites")




pdf("~/Desktop/heatmap.pdf")
give_heatmap(carb)
give_heatmap(ene)
give_heatmap(lipid)
give_heatmap(bile)
give_heatmap(nuc)
give_heatmap(aa)
give_heatmap(oth_aa)
give_heatmap(glycan)
give_heatmap(co_vit)
give_heatmap(terp)
give_heatmap(sec_met)
dev.off()

