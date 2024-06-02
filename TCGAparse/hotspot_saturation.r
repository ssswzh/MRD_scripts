#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# @Time    : 2022/01/07
# @Author  : zhangsiwen
# @contact : zhang.siwen@puruijizhun.com
# @LastModified: 2022/01/14
# @ChangeLog
#     20220107, first version


library(optparse)
library(ggplot2)
theme_set(theme_bw())


# passing arguments
option_list = list(
    make_option("--file", type="character", default=NULL, 
              help="dataset file name", metavar="character"),
    make_option("--out", type="character", default="file.saturation.png", 
              help="output file name [default= %default]", metavar="character")
); 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("Please assign input and output file.n", call.=FALSE)
}


# read dataframe
#ma <- read.table("F:\\MRD\\TCGA\\test.Not_Applicable.HotSites.tsv", sep="\t", header=TRUE, row.names=1)
ma <- read.table(opt$file, sep="\t", header=TRUE, row.names=1, comment.char = "#")
ma <- ma[, c("Covered_Sample_Number","Cumulate_Sample_Number")]
ma$Site <- 1:dim(ma)[1]
# total sample number
total_sample <- ma$Cumulate_Sample_Number[dim(ma)[1]]


# extract partial dataframe
if(dim(ma)[1]>=1000){
    site50 <- paste0(round(ma[50,'Cumulate_Sample_Number']/total_sample*100,2),"%")
    site100 <- paste0(round(ma[100,'Cumulate_Sample_Number']/total_sample*100,2),"%")
    site200 <- paste0(round(ma[200,'Cumulate_Sample_Number']/total_sample*100,2),"%")
    site300 <- paste0(round(ma[300,'Cumulate_Sample_Number']/total_sample*100,2),"%")
    site500 <- paste0(round(ma[500,'Cumulate_Sample_Number']/total_sample*100,2),"%")
    site800 <- paste0(round(ma[800,'Cumulate_Sample_Number']/total_sample*100,2),"%")
    site1000 <- paste0(round(ma[1000,'Cumulate_Sample_Number']/total_sample*100,2),"%")
} else if(dim(ma)[1]>=500) {
    site50 <- paste0(round(ma[50,'Cumulate_Sample_Number']/total_sample*100,2),"%")
    site100 <- paste0(round(ma[100,'Cumulate_Sample_Number']/total_sample*100,2),"%")
    site200 <- paste0(round(ma[200,'Cumulate_Sample_Number']/total_sample*100,2),"%")
    site300 <- paste0(round(ma[300,'Cumulate_Sample_Number']/total_sample*100,2),"%")
    site500 <- paste0(round(ma[500,'Cumulate_Sample_Number']/total_sample*100,2),"%")
    site800 <- ""
    site1000 <- ""
} else if(dim(ma)[1]>=200) {
    site50 <- paste0(round(ma[50,'Cumulate_Sample_Number']/total_sample*100,2),"%")
    site100 <- paste0(round(ma[100,'Cumulate_Sample_Number']/total_sample*100,2),"%")
    site200 <- paste0(round(ma[200,'Cumulate_Sample_Number']/total_sample*100,2),"%")
    site300 <- ""
    site500 <- ""
    site800 <- ""
    site1000 <- ""
}
ymax <- max(ma$Cumulate_Sample_Number)
ymin1 <- ymax/10
ymin2 <- ymax/20


# plot
ggplot(ma[1:1000,], aes(x=Site)) + 
    geom_line(aes(y=Covered_Sample_Number), color="#00ba38") + geom_text(aes(x=1000, y=0, label="Cover"), color="#00ba38", size=3) +
    geom_line(aes(y=Cumulate_Sample_Number), color="#f8766d") + geom_text(aes(x=1000, y=ymax-ymin2, label="Cumulate"), color="#f8766d", size=3) +
    labs(title="Saturation Curve for Detected Mutation Sites and Covered Sample Number", 
        #caption="test", y="Sample Number", size=3) + 
        caption=paste0("Source: ", opt$file), y="Sample Number", size=3) + # title and caption
    scale_x_continuous(breaks=seq(0,1000,100),limits=c(0,1000)) + 
    #geom_vline(xintercept=c(50,100,200,300,500,800,1000), linetype='dashed', color='blue', size=1) +
    geom_text(aes(x=50, y=ymin2, label=site50), color="black", size=3) + 
    geom_text(aes(x=100, y=ymin1, label=site100), color="black", size=3) + 
    geom_text(aes(x=200, y=ymin2, label=site200), color="black", size=3) + 
    geom_text(aes(x=300, y=ymin1, label=site300), color="black", size=3) + 
    geom_text(aes(x=500, y=ymin2, label=site500), color="black", size=3) + 
    geom_text(aes(x=800, y=ymin1, label=site800), color="black", size=3) + 
    geom_text(aes(x=1000, y=ymin2, label=site1000), color="black", size=3) + 
    theme(panel.grid.major = element_line(color = "grey", size = 0.5, linetype = 'dashed')) + 
    theme(plot.title=element_text(size=10), axis.text=element_text(size=10), axis.title=element_text(size=10)) + 
ggsave(filename=paste(opt$out,"saturation.png",sep="."), width=16, height=9, units='in', dpi=300)

# export saturation stats
stats <- data.frame("SitesCovered"=c(50,100,200,300,500,800,1000), "SamplesCovered"=c(site50,site100,site200,site300,site500,site800,site1000), row.names = NULL, stringsAsFactors = TRUE)
write.table(stats, paste(opt$out,"saturation.tsv",sep="."), quote = FALSE, sep="\t", row.names = FALSE)
