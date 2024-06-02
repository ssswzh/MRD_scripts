#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# @Time    : 2022/08/23
# @Author  : zhangsiwen
# @contact : zhang.siwen
# @LastModified: 2022/08/23
# @ChangeLog
#     20220823, first version
#     20220919, add --plot, change DrawScatterPlot() to draw multiple samples in a single file
#     20221229, add --overallvaf, change overall vaf to a optional argument

.libPaths(c("/mnt/ddngs/zhangsw/miniconda3/envs/R4.2.0/lib/R/library", "/mnt/ddngs/zhangsw/miniconda3/lib/R/library"))

# iDES polishing

library(fitdistrplus)
library(RColorBrewer)
library(data.table)
library(stringr)
library(tidyr)
library(BSDA)
library(optparse)
library(ggplot2)
theme_set(theme_classic())


# read file with format: mutationid, sample1, sample2,..., sampleN
ProcessData <- function ( inputfile ) {
  
  df <- read.table(inputfile, header=T, sep="\t")
  rownames(df) <- df$mutation_id
  nctcols <- colnames(df)[2:dim(df)[2]]
  df[,nctcols] <- as.numeric(unlist(df[,nctcols]))
  
  df <- separate(data=df, col=mutation_id, sep=":", into=c('chrom','pos','ref','alt'))
  df$NoiseType <- paste0(df$ref, '>', df$alt)
  df$id <- paste0(df$chrom, ':', df$pos)
  
  #df <- df[str_length(df$ref)==1 & str_length(df$alt)==1,]
  
  return(list(df, nctcols))
  
}


# calculate mean, std, fit data using weibull distribution by row
FitModelbyType <- function ( df, nctcols, scale=10000, remove_zero=T ) {
  
  df$delta <- apply(df[,nctcols] != 0, 1, sum)
  if( remove_zero==T ) {
    df <- df[df$delta!=0, ]
  }
  
  new_df <- copy(df)
  value_cols <- c("Mean","Median","STD","FitShape","FitScale","Correlation","CorPvalue")
  
  for ( i in rownames(df) ) {
    mean <- mean(as.numeric(df[i,nctcols]))
    median <- median(as.numeric(df[i,nctcols]))
    std <- sd(as.numeric(df[i,nctcols]))
    fit <- NA
    tryCatch({
      fit <- suppressWarnings(summary(fitdist(as.numeric(df[i,nctcols])*scale, "weibull"))$estimate)}, 
      warning = function(war){ b<-1 }, error = function(err){ fit <<- NA }) 
    #fit <- fitdist(as.numeric(df[i,nctcols]*scale), 'weibull')
    # correlation
    if ( class(fit) == "logical" ) {
      co <- data.frame(estimate=NA, p.value=NA)
    } else {
      w <- qqplot(as.numeric(df[i,nctcols])*scale, rweibull(100, fit['shape'], fit['scale']))
      co <- cor.test(w$x, w$y)      
    }

    new_df[i,value_cols] <- c(mean, median, std, fit['shape'], fit['scale']/scale, co$estimate, co$p.value)
  }
  
  return(new_df)
  
}


# get vaf cutoff by zscore
SolveVAFbyZscore <- function ( zscore, mean, std ) {
  return(zscore * std + mean)
}


# get vaf by weibull distribution
SolveVAFbyWeibull <- function ( pvalue, shape, scale, delta, corrected=T ) {
  if ( corrected==T ) {
    # pvalue <- 1-((1-delta)+(delta*p))
    pvalue <- (delta-pvalue) / delta
  }
  return(qweibull(pvalue, shape, scale, lower.tail=TRUE))
}


# get overall background vaf mean and std
OverallVAFbyGaussian <- function ( df, nctcols, remove_zero=F, pvalue=0.05, maxvaf=0.01 ){

  df$delta <- apply(df[,nctcols] != 0, 1, sum)
  if( remove_zero==T ) {
    df <- df[df$delta!=0, ]
  }
  
  # filter vaf 
  df$maxvaf <- apply(df[,nctcols], 1, max)
  df <- df[df$maxvaf<=maxvaf, ]
  df$maxvaf <- NULL
  
  bg <- dim(df)[1]
  zscore <- qnorm(1 - (pvalue/bg))
  all_data <- as.numeric(unlist(df[,nctcols]))
  #delta <- sum(all_data!=0)
  # fit <- suppressWarnings(
  #   summary(
  #     fitdist(all_data*scale, "weibull", lower=c(0,0), start=list(shape=mean(all_data*scale), scale=sd(all_data*scale)))
  #     )$estimate)
  
  mean <- mean(all_data)
  std <- sd(all_data)
  vaf_cutoff <- SolveVAFbyZscore(zscore, mean, std)

  return(vaf_cutoff)

}


# calculate VAF by gaussian or weibull
ModelCalculateVAFCutoff <- function ( model_df, pvalue=0.05, vaf_cutoff=0.0003 ) {
  
  # P-value cutoffs with Bonferroni correction
  bg <- dim(model_df)[1]
  zscore <- qnorm(1 - (pvalue/bg))
  weibullp <- pvalue/bg
  
  for ( i in rownames(model_df) ) {
    
    if ( is.na(model_df[i,'FitShape'] )) {
      # calculate vaf by zscore, mean and std
      model_df[i,'VAFcutoff'] <- SolveVAFbyZscore(zscore, model_df[i,'Mean'], model_df[i,'STD'])
      if ( model_df[i,'VAFcutoff']==0 ) {
        model_df[i,'VAFcutoff'] <- vaf_cutoff
      }
    }
    else {
      # calculate vaf by weibull
      model_df[i,'VAFcutoff'] <- SolveVAFbyWeibull(pvalue, model_df[i,'FitShape'], model_df[i,'FitScale'], model_df[i,'delta'], corrected=T) 
    }
    
  }
  
  return(model_df)
  
}


# test polishing result
ModelTest <- function ( model_df, sample_df, sample_cols, af_cutoff=0.0003 ) {
  
  new_df <- copy(sample_df)
  for ( sample in sample_cols ) {
    
    res_col <- paste("Result", sample, sep=':')
    for ( row in rownames(new_df) ) {
      
      if ( new_df[row,sample]==0 ) {
        new_df[row,res_col] <- 'Negative'
      } 
      else {
        if ( row %in% rownames(model_df) ) {
          if ( model_df[row,'VAFcutoff']==0 ) {
            new_df[row,res_col] <- ifelse( new_df[row,sample] >= af_cutoff, 'Positive', 'Negative')
          }
          else {
            new_df[row,res_col] <- ifelse( new_df[row,sample] >= model_df[row,'VAFcutoff'], 'Positive', 'Negative')
          }
        }
        else {
          new_df[row,res_col] <- ifelse( new_df[row,sample] >= af_cutoff, 'Positive', 'Negative')
        }
      }

    }
    
  }
  return(new_df)
  
}


# merge depth and vaf for ploting scatter plot
ExpandMergeDPVAF <- function ( depthDF, vafDF, col ) {
  
  depthDF <- depthDF[rownames(depthDF) %in% rownames(vafDF), ]
  
  # depth and freq matrix
  expandDf <- data.frame(matrix(nrow=dim(depthDF)[1]*length(col), ncol=5))
  colnames(expandDf) <- c('mutation_id','Depth', 'Frequency', 'Sample', 'NoiseType')
  for (i in 1:length(col)) {
    expandDf[((i-1)*dim(depthDF)[1]+1):(i*dim(depthDF)[1]),'mutation_id'] <- rownames(depthDF)
    expandDf[((i-1)*dim(depthDF)[1]+1):(i*dim(depthDF)[1]),'Depth'] <- depthDF[,col[i]]
    expandDf[((i-1)*dim(depthDF)[1]+1):(i*dim(depthDF)[1]),'Frequency'] <- vafDF[,col[i]]
    expandDf[((i-1)*dim(depthDF)[1]+1):(i*dim(depthDF)[1]),'Sample'] <- col[i]
    expandDf[((i-1)*dim(depthDF)[1]+1):(i*dim(depthDF)[1]),'NoiseType'] <- depthDF[,'NoiseType']
  }
  return(expandDf)
  
}


# main function
Main <- function (nct_file, test_file, outfile, overallvaf=FALSE, pvalue=0.05) {
  
  # build model
  
  NCTdata <- ProcessData(nct_file)
  NCTdf <- NCTdata[[1]]
  NCTcols <- NCTdata[[2]]
  #af_vector <- seq(0.00001, 0.01, 0.00001)
  
  NCTfit <- FitModelbyType(NCTdf, NCTcols, scale=10000, remove_zero=F)
  if ( overallvaf == TRUE) { # 20221229
    overall_bgaf <- OverallVAFbyGaussian(NCTdf, NCTcols, remove_zero=T, pvalue=numeric(pvalue), maxvaf=0.01)
  } else {
    overall_bgaf <- 0
  }
  NCTmodel <- ModelCalculateVAFCutoff(NCTfit, pvalue=numeric(pvalue), vaf_cutoff=overall_bgaf)
  
  # test model
  
  TESTdata <- ProcessData(test_file)
  TESTdf <- TESTdata[[1]]
  TESTcols <- TESTdata[[2]]
  TESTresult <- ModelTest(NCTmodel, TESTdf, TESTcols, af_cutoff=overall_bgaf)
  
  #outfile <- paste(sub('\\.tsv$', '', test_file), 'polishing.tsv', sep='.')
  write.table(NCTmodel, paste0(outfile,'.model_fit.tsv'), quote=F, sep="\t", row.names=F)
  write.table(TESTresult, paste0(outfile,'.polishing.tsv'), quote=F, sep="\t", row.names=F)
  
  print(overall_bgaf)
  
  dfs <- list(NCTdf, NCTcols, NCTmodel, TESTdf, TESTcols, TESTresult, overall_bgaf)
  names(dfs) <- c("NCTdf", "NCTcols", "NCTmodel", "TESTdf", "TESTcols", "TESTresult", "overall_bgaf")
  return(dfs)
  
}


# sites scatter plot
DrawScatterPlot <- function( NCTvaf, TESTvaf, NCTdp, TESTdp, out, plot=NULL, ymax=0.005 ) {
  
  NCTdepthdata <- ProcessData(NCTdp)
  NCTdepthdf <- NCTdepthdata[[1]]
  NCTdepthcol <- NCTdepthdata[[2]]
  
  TESTdepthdata <- ProcessData(TESTdp)
  TESTdepthdf <- TESTdepthdata[[1]]
  TESTdepthcol <- TESTdepthdata[[2]]
  
  # depth and vaf
  NCTscatterDF <- ExpandMergeDPVAF(NCTdepthdf, NCTvaf, NCTdepthcol)
  TESTscatterDF <- ExpandMergeDPVAF(TESTdepthdf, TESTvaf, TESTdepthcol)
  for (i in unique(TESTscatterDF$mutation_id) ) {
    for ( j in unique(TESTscatterDF$Sample) ) {
      TESTscatterDF[TESTscatterDF$mutation_id==i&TESTscatterDF$Sample==j,'Polishing'] <- TESTvaf[i, paste0('Result:',j)]
    }
  }
  
  # for values larger than ymax
  TESTscatterDF[TESTscatterDF$Frequency>=ymax, 'Frequency'] <- ymax
  
  # ad curve
  adcurve <- data.frame(Depth=seq(0,5e5,100), AD=10)
  adcurve$freq <- adcurve$AD/adcurve$Depth
  
  # scatter plot 20220919
  if (plot == "all") {
    pdf(paste0(out, '.scatterplot.pdf'), width=8, height=6)
    print(
      ggplot(NCTscatterDF, aes(x=Depth, y=Frequency), color='black') + 
        geom_point() +
        geom_hline(yintercept=0.0003, color='black') + 
        geom_text(x=4.5e5, y=0.0003, label='y=0.0003', color='black', vjust=-0.5) +
        geom_hline(yintercept=0.0005, color='black') + 
        geom_text(x=4.5e5, y=0.0005, label='y=0.0005', color='black', vjust=-0.5) +
        geom_hline(yintercept=0.0001, color='black') + 
        geom_text(x=4.5e5, y=0.0001, label='y=0.0001', color='black', vjust=-0.5) +
        geom_vline(xintercept=1000, color='black') + 
        geom_text(x=0, y=0.005, label='x=1000', color='black', hjust=-0.5, vjust=0.5) +
        ylim(0,ymax) +
        geom_line(data=adcurve, aes(x=Depth, y=freq), color='black', linewidth=1) +
        geom_text(x=4.5e5, y=0, label='AD=10', color='black', vjust=1) +
        geom_point(data=TESTscatterDF, aes(x=Depth, y=Frequency, color=Polishing, shape=Sample), size=3) +
        scale_color_brewer(palette='Set1')
    )
    dev.off()
  } else {
    for (sam in TESTdepthcol) {
    pdf(paste(out, sam, 'scatterplot.pdf', sep='.'), width=8, height=6)
    print(
      ggplot(NCTscatterDF, aes(x=Depth, y=Frequency), color='black') + 
        geom_point() +
        geom_hline(yintercept=0.0003, color='black') + 
        geom_text(x=4.5e5, y=0.0003, label='y=0.0003', color='black', vjust=-0.5) +
        geom_hline(yintercept=0.0005, color='black') + 
        geom_text(x=4.5e5, y=0.0005, label='y=0.0005', color='black', vjust=-0.5) +
        geom_hline(yintercept=0.0001, color='black') + 
        geom_text(x=4.5e5, y=0.0001, label='y=0.0001', color='black', vjust=-0.5) +
        geom_vline(xintercept=1000, color='black') + 
        geom_text(x=0, y=0.005, label='x=1000', color='black', hjust=-0.5, vjust=0.5) +
        ylim(0,ymax) +
        geom_line(data=adcurve, aes(x=Depth, y=freq), color='black', size=1) +
        geom_text(x=4.5e5, y=0, label='AD=10', color='black', vjust=1) +
        geom_point(data=TESTscatterDF[TESTscatterDF$Sample==sam,], aes(x=Depth, y=Frequency, color=Polishing, shape=Polishing), size=3) +
        scale_color_brewer(palette='Set1')
    )
    dev.off()
    }
  }
  
  dfs <- list(NCTscatterDF, TESTscatterDF)
  names(dfs) <- c("NCTscatterDF", "TESTscatterDF")
  return(dfs)
  
}



# main run

option_list = list(
  make_option("--nctvaf", type="character", default=NULL, help="nct vaf tsv file, first column is mutation_id, other columns are VAF of different NCT samples"),
  make_option("--nctdp", type="character", default=NULL, help="nct dp tsv file, first column is mutation_id, other columns are depth of different NCT samples"),
  make_option("--testvaf", type="character", default=NULL, help="test vaf tsv file, same format as --nctvaf file"),
  make_option("--testdp", type="character", default=NULL, help="test dp tsv file, same format as --nctdp file"),
  make_option("--plot", type="character", default=NULL, help="plot result by each sample or all sample, choose from [each, all],\nif do not want plot, don't specify this option"),
  make_option("--ylim", type="double", default=0.005, help="y-axis limit value for plot, default 0.005"),
  make_option("--overallvaf", type="character", default=FALSE, help="not calculate overall vaf for sites which have frequency=0 in all control samples, default FALSE"),
  make_option("--pval", type="double", default=0.05, help="not calculate overall vaf for sites which have frequency=0 in all control samples, default FALSE"),
  make_option("--out", type="character", default=NULL, help="output prefix")
);

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# list(NCTdf, NCTcols, NCTmodel, TESTdf, TESTcols, TESTresult, overall_bgaf)
results <- Main(opt$nctvaf, opt$testvaf, opt$out, opt$overallvaf)
if (is.null(opt$plot)) {
  print("Choose to not plot mutations .", call.=FALSE)
} else {
  dfs <- DrawScatterPlot(results$NCTdf, results$TESTresult, opt$nctdp, opt$testdp, opt$out, opt$plot, ymax=as.numeric(opt$ylim))
}

# results <- Main("nct.all_sample.VAF.tsv", "test.all_sample.VAF.tsv", "test.all_sample")
# dfs <- DrawScatterPlot(results$NCTdf, results$TESTresult, "nct.all_sample.DP.tsv", "test.all_sample.DP.tsv", "test.all_sample.polishing", 'all', ymax=0.005)

