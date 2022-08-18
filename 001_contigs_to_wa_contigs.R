# Install packages:
# install.packages("base", repos = "http://cran.us.r-project.org")
# install.packages("data.table", repos = "http://cran.us.r-project.org", dependencies = TRUE)
# install.packages("dplyr", repos = "http://cran.us.r-project.org")
# install.packages("seqinr", repos = "http://cran.us.r-project.org")
# install.packages("stringr", repos = "http://cran.us.r-project.org")
# install.packages("utils", repos = "http://cran.us.r-project.org")

#upload all libraries
library(base)
library(data.table)
library(dplyr)
library(seqinr)
library(stringr)
library(utils)

# input dir
pig.id.basedir = "/shared/homes/152324/out_new" #make sure you have permission to write in the subdir /metabat
## test here: 
#pig.id.basedir = "/Users/dgaio/cloudstor/Gaio/out_new_test"


# 1. gather all contig names (>k141") from fasta files and store them in a file "bins_to_contigs.csv"
for (pig.id in list.files(pig.id.basedir)) {
  pig.id.dir = file.path(pig.id.basedir, pig.id, "metabat")

  df = data.frame(
    pig = character(),
    bin = character(),
    contig = character()
  )
  
  fasta.files = list.files(pig.id.dir, pattern=".fa")
  for (fasta.file in fasta.files){
    ver <- read.fasta(
      file.path(pig.id.dir, fasta.file),
      as.string = TRUE,
      strip.desc = TRUE
    )
    annotation <- getAnnot(ver)
    #unlist to create a vector from a list (goes faster)
    annotation.vector = unlist(annotation)
    
    df = rbind(
      df,
      data.frame(
        pig = rep(pig.id, length(annotation.vector)),
        bin = rep(fasta.file, length(annotation.vector)),
        contig = annotation.vector
      )
    )
  }
  
  fwrite(
    x = df,
    file = file.path(pig.id.dir, "bins_to_contigs.csv"),
    row.names=FALSE
  )
}

# 2. this loop creates a file inside each sample directory from depth.txt files, stripping off 
#the sample IDs from all headers, leaving only sampling dates and extensions
for (pig.id in list.files(pig.id.basedir)) {
  pig.id.dir = file.path(pig.id.basedir, pig.id, "metabat")
  
  depth.files = list.files(pig.id.dir, pattern=".txt")
  for (depth.file in depth.files) {
    nub <- read.table(
      file.path(pig.id.dir, depth.file), 
      header=TRUE
    )
    colnames(nub)
    #function to strip the sample ID from the headers
    colClean <- function(nub){ colnames(nub) <- gsub("^[0-9a-zA-z]+_[0-9a-zA-z]+.", " ", colnames(nub)); nub }
    new <- colClean(nub)
  }
  fwrite(
    x = new,
    file = file.path(
      pig.id.dir, 
      "depth_clean.csv"
    ),
    row.names=FALSE
  )
}

# 3. merging depth_clean.csv with bins_to_contigs.csv
for (pig.id in list.files(pig.id.basedir)) {
  pig.id.dir = file.path(pig.id.basedir, pig.id, "metabat")
  
  infile1 <- paste(pig.id.dir, "depth_clean.csv", sep="/")
  infile2 <- paste(pig.id.dir, "bins_to_contigs.csv", sep="/")
  
  a <- read.table(
    infile1, 
    header=TRUE, 
    sep=",", 
    row.names=NULL
  )
  b <- read.table(
    infile2, 
    header=TRUE, 
    sep=",", 
    row.names=NULL
  )
  
  a_b <- merge(
    a, 
    b, 
    by.x="contigName", 
    by.y="contig"
  )
  
  fwrite(
    x = a_b,
    file = file.path(pig.id.dir, "depth_with_bins.csv"),
    row.names=FALSE
  )
}


#4. calculating weighted averages for contigs. Output: length-normalized contig depths
for (pig.id in list.files(pig.id.basedir)) {
  pig.id.dir = file.path(pig.id.basedir, pig.id, "metabat")
  
  old <- paste(pig.id.dir, "depth_with_bins.csv", sep="/")
  
  df <- read.table(
    file = old, 
    header = TRUE, 
    sep = ",", 
    row.names = NULL
  )
  
  df_new <- df %>% 
    select(-ends_with(".var"))
  dt <- as.data.table(df_new)
  colsTOweight <- dt %>% 
    select(contains('.bam'))
  colsTOweight <- colnames(colsTOweight)
  
  dt_new <- setDT(dt)[, lapply(.SD, `/`, dt$contigLen), .SDcols = colsTOweight]
  
  dt_new2 <- cbind(df$pig, df$bin,df$contigName,df$contigLen, dt_new)
  
  colnames(dt_new2)[1:4] <- c("pig","bin", "contigName","contigLen")
  
  fwrite(
    x = dt_new2,
    file = file.path(pig.id.dir, "wa_contigs.csv"),
    row.names=FALSE
  )
  
}
