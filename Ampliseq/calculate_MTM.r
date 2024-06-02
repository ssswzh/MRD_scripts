#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# @Author: zhang.siwen

library(optparse)

# Argument parsing
option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="polishing/all_sample.final.result.tsv"),
  make_option(c("-o", "--out"), type="character", default=NULL, 
              help="outpath/output.file"),
  make_option(c("-s", "--sample"), type="character", default=NULL,
              help="sample id used in analysis"),
  make_option(c("-m", "--mass"), type="double", default=NULL,
              help="cfDNAmass (ng)"),
  make_option(c("-p", "--plasma"), type="double", default=3,
              help="plasma volumn (mL), default '3'")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if ( FALSE ) {
  datafile <- "all_sample.final.result.tsv"
  id <- "T990500066BC01-112"
  cfDNAmass <- 29.05
  plasmaVolumn <- 3
}

calculateMeanVAF <- function ( df, id ) {
  id <- make.names(id)
  res <- paste0("Result.", id)
  tmp <- df[, c(id,res)]
  # negative site assign to 0
  tmp[tmp[,res]=="Negative", id] <- 0
  return ( mean(tmp[,id]) )
}

calculateMTM <- function (meanVAF, cfDNAmass, plasmaVolumn ) {
  return ( meanVAF * cfDNAmass * 1000 / (3.3 * plasmaVolumn) )
}



df <- read.table(opt$input, sep="\t", header=T)

MTM <- calculateMTM(calculateMeanVAF(df, opt$sample), opt$mass, opt$plasma)

write(MTM, opt$out)
