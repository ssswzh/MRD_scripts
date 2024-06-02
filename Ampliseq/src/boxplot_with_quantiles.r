#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# @Time    : 2022/07/25
# @Author  : zhangsiwen
# @contact : zhang.siwen
# @LastModified: 2022/09/19
# @ChangeLog
#     20220725, first version
#     20220919, allow multiple columns or single column to plot in a single boxplot


.libPaths(c("/mnt/ddngs/zhangsw/miniconda3/envs/R4.2.0/lib/R/library", "/mnt/ddngs/zhangsw/miniconda3/lib/R/library"))

usage <- "
  Usage:
    Rscript --vanilla boxplot_with_quantiles.r --input per_target_cov --out outprefix --dir ./ --key mean_coverage
    Rscript --vanilla boxplot_with_quantiles.r --input per_target_cov --out outprefix --dir ./
  
  For single per_target_cov file:
    must have keyword you specify by --key in first line, like output of Picard CollectTargetedPcrMetrics:
      e.g. (tab-delimited) (say --key=mean_coverage in this case) (other columns are ignored)
    chrom	start	end	length	name	%gc	mean_coverage	normalized_coverage	min_normalized_coverage	max_normalized_coverage	min_coverage max_coverage	pct_0x	read_count
    1	5965383	5965442	60	NPHP4(0.60625)	0.433333	123661.066667	0.276457	0.273001	0.277025	122115	123915	0  123919
    1	159558307	159558369	63	APCS(0.585429314830876)	0.492063	162237.746032	0.362699	0.35756	0.363319	15993  162515	0	162527
  
  For batch per_target_cov files:
    must have an identifier in first column, and other fields are sample values
      e.g. (tab-delimited) (do not specify --key) (or you can specify --key=NA12878 and will only plot NA12878 column)
    mutation_id	20220623012	NA12878	NA24143	NA24149	NA24694	NA24695
    1:40363293:C:G	60010	76252	72770	80851	78831	77868
    1:53250678:A:G	168148	183812	186746	206551	181213	207821
  
  Other values could be used as well, for instance, use VAF values instead of depth in the file.
  
"

# test if there is at least one argument: if not, return an error
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop(usage, call.=FALSE)
} else if (args[1]=='-h'||args[1]=='--help') {
  stop(usage, call.=FALSE)
}


library(optparse, quietly=TRUE)
library(ggplot2, quietly=TRUE)
#library(ggpubr, quietly=TRUE)
theme_set(theme_classic())
library(dplyr, quietly=TRUE)
library(data.table, quietly=TRUE)
library(RColorBrewer, quietly=TRUE)
library(ggpmisc, quietly=TRUE)


if (FALSE) {
  df <- read.table('LJ2103907BC01-7.per_target_cov', sep = '\t', header = T)
  out <- 'LJ2103907BC01-7.per_target_cov'
  dir <- './'
  key <- 'mean_coverage'
}


# coverage boxplot

SummarizeStats <- function ( df, id, vars ) {
  
  mat <- melt(df, id.vars=id, measure.vars=vars)
  
  mins <- mat %>% group_by(variable) %>% summarise_at(vars(value), list(Min = min))
  means <- mat %>% group_by(variable) %>% summarise_at(vars(value), list(Mean = mean))
  medians <- mat %>% group_by(variable) %>% summarise_at(vars(value), list(Median = median))
  maxs <- mat %>% group_by(variable) %>% summarise_at(vars(value), list(Max = max))
  stats <- data.frame(mins, medians, means, maxs)[,c(1,2,4,6,8)]
  colnames(stats) <- c('variable', 'Min', 'Median','Mean', 'Max')
  
  if (max(stats$Max)>=1) {
    stats <- stats %>% mutate_if(is.numeric, round)
  } else {
    stats <- stats %>% mutate_if(is.numeric, round, 4)
  }
  
  return(stats)
  
}


DrawBoxplot <- function (df, stats, key, out, width=4, height=6) {
  
  pdf(paste0(out,'.boxplot.pdf'), width = width, height = height)
  print(
    ggplot(df, aes_string(x='variable', y=key, color='variable', fill='variable')) +
      geom_boxplot(alpha=0.5, outlier.color=NA) + geom_jitter(shape=16, position=position_jitter(0.2), size=2) + 
      geom_text(data = stats, aes(x = variable, y = Min, label = paste('Min:',Min)), size = 4, color='black') + 
      geom_text(data = stats, aes(x = variable, y = Median, label = paste('Median:',Median)), size = 4, color='black') + 
      geom_text(data = stats, aes(x = variable, y = Mean, label = paste('Mean:',Mean)), size = 4, color='black') + 
      geom_text(data = stats, aes(x = variable, y = Max, label = paste('Max:',Max)), size = 4, color='black') + 
      labs(x='Groups', y=key) + 
      #scale_y_continuous(trans='log10') +
      #scale_color_brewer(palette="Set2") + scale_fill_brewer(palette="Set2") +
      theme(axis.text=element_text(size=10), axis.title=element_text(size=10), 
            plot.title=element_text(size=15), plot.caption=element_text(size=10), 
            #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
            legend.position="none")
  )
  dev.off()
  
}


Main <- function (inputfile, outdir='./', outfile, key=NULL) {
  
  df <- read.table(inputfile, sep = '\t', header = T)
  if (outdir != './') {
    dir.create(file.path(outdir), showWarnings = FALSE)
    outfile <- paste(outdir, outfile, sep='/')
  }
  
  if (is.null(key)) { # multiple column values
    width <- 2+dim(df)[2]
    height <- 6
    sample_name <- colnames(df)[2:dim(df)[2]]
    stats <- SummarizeStats(df, colnames(df)[1], vars=sample_name)
    df <- melt(df, id.vars=colnames(df)[1], measure.vars=sample_name)
    DrawBoxplot(df, stats, 'value', outfile, width=width, height=height)
    
  } else { # single value with key name
    width <- 3
    height <- 6
    sample_name <- strsplit(basename(inputfile), '[.]')[[1]][1]
    df$Sample <- sample_name
    df[,key] <- as.numeric(df[,key])
    #df <- df[,c('Sample',key)]
    stats <- SummarizeStats(df, id='Sample', vars=key) 
    stats$variable <- sample_name
    df$variable <- sample_name
    DrawBoxplot(df, stats, key, outfile, width=width, height=height)
  }
  
  return(stats)
  
}



# define arguments
option_list = list(
  make_option("--input", type="character", default=NULL, help="input file, details see usage", metavar="character"),
  make_option("--out", type="character", default=NULL, help="output prefix", metavar="character"),
  make_option("--dir", type="character", default='./', help="output outdir, default current path", metavar="character"),
  make_option("--key", type="character", default=NULL, help="", metavar="character")
);

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

stats <- Main(opt$input, opt$dir, opt$out, opt$key)
# stats <- Main("LJ2103907BC01-7.per_target_cov", outdir='./', "LJ2103907BC01-7.per_target_cov", key='mean_coverage')
# stats <- Main("../polishing/SAM001_kit/polishing/nct.all_sample.DP.tsv", outdir='./', "../polishing/SAM001_kit/polishing/nct.all_sample.DP")

