#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# @Time    : 2022/01/17
# @Author  : zhangsiwen
# @contact : zhang.siwen@puruijizhun.com
# @LastModified: 2022/01/17
# @ChangeLog
#     20220117, first version, modified from hotspot_saturation.r


library(optparse)
library(ggplot2)
library(RColorBrewer)
theme_set(theme_bw())


# passing arguments
option_list = list(
    make_option("--file", type="character", default=NULL, 
        help="project name and dataset file separated by Tab, first line will be Project\tFile", metavar="character"),
    make_option("--out", type="character", 
        help="output file name, will add saturation.png", default=NULL, metavar="character"),
    make_option("--num", type="character", 
        help="number of sites to be plot, default 1000", default="1000", metavar="character")
); 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("Please assign input and output file.n", call.=FALSE)
}



# read dataframe
#files <- read.table("C:\\Users\\zhangsw\\Desktop\\test.tsv", header=TRUE, sep="\t")
files <- read.table(opt$file, header=TRUE, sep="\t")
# num <- 1000
num <- as.numeric(opt$num)

for(i in 1:dim(files)[1]) {
    ma <- read.table(files$File[i], header=TRUE, sep="\t")
    if(files$Project[i]=="TCGA") {
        ma$Cumulate_Sample_Number_Percent <- ma$Cumulate_Sample_Number/ma$Cumulate_Sample_Number[dim(ma)[1]]
    }
    if("Cumulate_Sample_Number_Percent" %in% colnames(ma)) {
        tmp <- data.frame("Sites"=1:num, "Sample"=ma$Cumulate_Sample_Number_Percent[1:num], "Project"=rep(files$Project[i],num), row.names = NULL, stringsAsFactors = TRUE)
    }
    if("Dataset_Cumulate_Sample_Number_Percent" %in% colnames(ma)) {
        tmp <- data.frame("Sites"=1:num, "Sample"=ma$Dataset_Cumulate_Sample_Number_Percent[1:num], "Project"=rep(files$Project[i],num), row.names = NULL, stringsAsFactors = TRUE)
    }
    if(i==1) {
        mergedDf <- tmp
    } else {
        mergedDf <- rbind(mergedDf, tmp)
    }
}
mergedDf$Percent <- round(mergedDf$Sample*100,2)
#mid <- mergedDf[which(mergedDf$Sites==num/2),'Percent']

ggplot(mergedDf, aes(x=Sites, y=Percent, color=Project, fill=Project)) + geom_line(size=1) + 
    labs(title="Mutation Covered Sample Percent In Datasets", 
        y="Covered Sample Percent (%)", size=3) + # title and caption
    scale_y_continuous(breaks=seq(0,100,10),limits=c(0,100)) +
    scale_x_continuous(breaks=seq(0,num,100),limits=c(0,num)) +
    #scale_color_jama() + scale_fill_jama() + 
    scale_color_brewer(palette = "Set2") + 
    theme(panel.grid.major = element_line(color = "grey", size = 0.5, linetype = 'dashed')) + 
    theme(plot.title=element_text(size=10), axis.text=element_text(size=10), axis.title=element_text(size=10)) + 
ggsave(filename=paste(opt$out,"saturation.png",sep="."), width=8, height=4.5, units='in', dpi=300)

# output stats
s <- c(1,2,3,4,5,6,7,8,9,10,20,30,50,100,200,300,500,800,1000)
stats <- mergedDf[which(mergedDf$Sites %in% s),c('Sites','Percent','Project')]
write.table(stats, paste(opt$out,"saturation.tsv",sep="."), quote = FALSE, sep="\t", row.names = FALSE)
