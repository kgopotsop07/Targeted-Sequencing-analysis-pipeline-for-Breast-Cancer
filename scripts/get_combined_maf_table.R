#!/usr/bin/env Rscript

library(maftools)
library(dplyr)

# get input argument
args = commandArgs(trailingOnly=TRUE)

#stop the script if no command line argument

if(length(args)==0){
  print("Please include differential expression results!")
  stop("Requires command line argument.")
  } else {
  ## specifying the input arguments
  path.in <- args[1]
  reference_type <- args[2]
}
##----------------------
# install and load libraries
#r = getOption("repos")
#r["CRAN"] = "http://cran.us.r-project.org"
#options(repos = r)

#packages = c("dada2", "Biostrings",
#             "ShortRead", "ggplot2", "reshape2")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

## Now load or install&load all
#package.check <- lapply(
#  packages,
#  FUN = function(x) {
#    if (!require(x, character.only = TRUE)) {
#      BiocManager::install(x, dependencies = TRUE)
#      library(x, character.only = TRUE)
#    }
#  }
#)

#invisible(package.check)
## filter to use non empty files
setwd(as.character(path.in))

variant_annotation_list <-sort (list.files(pattern = "_multianno.txt"))
info <- file.info(variant_annotation_list)
file_list <- rownames(info[info$size != 0L, ])
TSVList<-lapply(file_list ,function(x){read.table(x,header=T,sep = "\t")})
TSVs<-gsub(file_list,pattern = "_multianno.txt",replacement = "")
names(TSVList)<-TSVs

## input files as mafs
MAFs<-sort(list.files(pattern = ".hg38_multianno.txt"))
MAFs_info <- file.info(variant_annotation_list)
MAFs_file_list <- rownames(info[info$size != 0L, ])

MafList<-lapply(MAFs_file_list,function(x){annovarToMaf(x,Center = NULL,
   refBuild = as.character(reference_type) ,tsbCol = NULL,
   table = "refGene",ens2hugo = TRUE,
   basename = NULL,sep = "\t",
   MAFobj = FALSE,sampleAnno = NULL)})

names(MafList)<-TSVs

## writing output files
docs <- bind_rows(MafList)
write.table(docs,paste0("multiannotation_joined_output.tsv"), sep = "\t")
