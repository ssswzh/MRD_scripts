#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# @Time    : 2022/01/14
# @Author  : zhangsiwen
# @contact : zhang.siwen@puruijizhun.com
# @LastModified: 2022/01/14
# @ChangeLog
#     20220114, first version, modified from hotspot_saturation.r


library(optparse)
library(ggplot2)
theme_set(theme_bw())


# passing arguments
option_list = list(
    make_option("--file", type="character", default=NULL, 
        help="dataset file name", metavar="character"),
    make_option("--out", type="character", 
        help="output file name, will add saturation.png", metavar="character"),
    make_option("--baseline", type="character", default="Baseline_Cumulate_Sample_Number_Percent", 
        help="baseline column name [default=%default]", metavar="character"),
    make_option("--dataset", type="character", default="Dataset_Cumulate_Sample_Number_Percent", 
        help="dataset column name [default=%default]", metavar="character")
); 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("Please assign input and output file.n", call.=FALSE)
}


# read dataframe
#ma <- read.table("F:\\MRD\\cBioPortal\\data_mutation.Compare.BRCA.HotSites.tsv", sep="\t", header=TRUE, comment.char = "#")
ma <- read.table(opt$file, sep="\t", header=TRUE, comment.char = "#")
baseline <- opt$baseline
dataset <- opt$dataset
ma$Baseline <- ma[,baseline] * 100
ma$Dataset <- ma[,dataset] * 100
ma$Site <- 1:dim(ma)[1]

base_site50 <- round(ma[50,'Baseline'],2)
base_site100 <- round(ma[100,'Baseline'],2)
base_site200 <- round(ma[200,'Baseline'],2)
base_site300 <- round(ma[300,'Baseline'],2)
base_site500 <- round(ma[500,'Baseline'],2)
base_site800 <- round(ma[800,'Baseline'],2)
base_site1000 <- round(ma[1000,'Baseline'],2)

data_site50 <- round(ma[50,'Dataset'],2)
data_site100 <- round(ma[100,'Dataset'],2)
data_site200 <- round(ma[200,'Dataset'],2)
data_site300 <- round(ma[300,'Dataset'],2)
data_site500 <- round(ma[500,'Dataset'],2)
data_site800 <- round(ma[800,'Dataset'],2)
data_site1000 <- round(ma[1000,'Dataset'],2)


# plot
ggplot(ma[1:1000,], aes(x=Site)) + 
    geom_line(aes(y=Baseline), color="#00ba38") + geom_text(aes(x=1000, y=10, label=baseline), color="#00ba38", size=3) +
    geom_line(aes(y=Dataset), color="#f8766d") + geom_text(aes(x=1000, y=5, label=dataset), color="#f8766d", size=3) +
    labs(title="Mutation Covered Sample Percent In Baseline and Dataset", 
        #caption="test", y="Covered Sample Percent (%)", size=3) + 
        caption=paste0("Source: ", opt$file), y="Covered Sample Percent (%)", size=3) + # title and caption
    scale_x_continuous(breaks=seq(0,1000,100),limits=c(0,1000)) + 
    #geom_vline(xintercept=c(50,100,200,300,500,800,1000), linetype='dashed', color='blue', size=1) +
    geom_text(aes(x=50, y=base_site50+5, label=paste0(base_site50,"%")), color="#00ba38", size=3) + 
    geom_text(aes(x=100, y=base_site100+5, label=paste0(base_site100,"%")), color="#00ba38", size=3) + 
    geom_text(aes(x=200, y=base_site200+5, label=paste0(base_site200,"%")), color="#00ba38", size=3) + 
    geom_text(aes(x=300, y=base_site300+5, label=paste0(base_site300,"%")), color="#00ba38", size=3) + 
    geom_text(aes(x=500, y=base_site500+5, label=paste0(base_site500,"%")), color="#00ba38", size=3) + 
    geom_text(aes(x=800, y=base_site800+5, label=paste0(base_site800,"%")), color="#00ba38", size=3) + 
    geom_text(aes(x=1000, y=base_site1000+5, label=paste0(base_site1000,"%")), color="#00ba38", size=3) + 
    geom_text(aes(x=50, y=data_site50-5, label=paste0(data_site50,"%")), color="#f8766d", size=3) + 
    geom_text(aes(x=100, y=data_site100-5, label=paste0(data_site100,"%")), color="#f8766d", size=3) + 
    geom_text(aes(x=200, y=data_site200-5, label=paste0(data_site200,"%")), color="#f8766d", size=3) + 
    geom_text(aes(x=300, y=data_site300-5, label=paste0(data_site300,"%")), color="#f8766d", size=3) + 
    geom_text(aes(x=500, y=data_site500-5, label=paste0(data_site500,"%")), color="#f8766d", size=3) + 
    geom_text(aes(x=800, y=data_site800-5, label=paste0(data_site800,"%")), color="#f8766d", size=3) + 
    geom_text(aes(x=1000, y=data_site1000-5, label=paste0(data_site1000,"%")), color="#f8766d", size=3) + 
    theme(panel.grid.major = element_line(color = "grey", size = 0.5, linetype = 'dashed')) + 
    theme(plot.title=element_text(size=10), axis.text=element_text(size=10), axis.title=element_text(size=10)) + 
ggsave(filename=paste(opt$out,"saturation.png",sep="."), width=8, height=4.5, units='in', dpi=300)


# export saturation stats
stats <- data.frame("SitesCovered"=c(50,100,200,300,500,800,1000), "BaselineSamplesCovered(%)"=c(base_site50,base_site100,base_site200,base_site300,base_site500,base_site800,base_site1000), "DatasetSamplesCovered(%)"=c(data_site50,data_site100,data_site200,data_site300,data_site500,data_site800,data_site1000), row.names = NULL, stringsAsFactors = TRUE)
write.table(stats, paste(opt$out,"saturation.tsv",sep="."), quote = FALSE, sep="\t", row.names = FALSE)
