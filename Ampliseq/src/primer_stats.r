#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0 || args[1]=="-h" || args[1]=="--help") {
  stop("\nUsage:
  Rscript --vanilla primer_stats.r sorted.bam.mapping_position.count.target.max_count.primer.tsv
  ", call.=FALSE)
}


library(ggplot2, quietly=TRUE)
theme_set(theme_classic())
library(RColorBrewer, quietly=TRUE)
library(dplyr, quietly=TRUE)


if (FALSE) {
  filename <- "P11.sorted.bam.mapping_position.count.target.max_count.primer.tsv"
  df <- read.table(filename, sep='\t',header=T)
  dir <- dirname(filename)
  base <- basename(sub('\\.tsv$', '', filename))
  out <- paste(dir, base, sep='/')
}


df <- read.table(args[1], sep='\t',header=T)
dir <- dirname(args[1])
base <- basename(sub('\\.tsv$', '', args[1]))
out <- paste(dir, base, sep='/')


# size dist

mins <- df %>% group_by(primer_type) %>% summarise_at(vars(seq_len), list(name = min))
means <- df %>% group_by(primer_type) %>% summarise_at(vars(seq_len), list(name = mean))
medians <- df %>% group_by(primer_type) %>% summarise_at(vars(seq_len), list(name = median))
maxs <- df %>% group_by(primer_type) %>% summarise_at(vars(seq_len), list(name = max))
stats <- data.frame(mins, medians, means, maxs)[,c(1,2,4,6,8)]
colnames(stats) <- c('primer_type', 'Min', 'Median','Mean', 'Max')
stats <- stats %>% mutate_if(is.numeric, round)

pdf(paste0(out,'.primer_size.pdf'), width = 6, height = 6)
ggplot(df, aes(x=primer_type, y=seq_len, color=primer_type, fill=primer_type)) +
  geom_boxplot(alpha=0.5, outlier.color=NA) + geom_jitter(shape=16, width=0.2, height=0, size=2) + 
  geom_text(data = stats, aes(x = primer_type, y = Min, label = paste('Min:',Min)), size = 4, color='black') + 
  geom_text(data = stats, aes(x = primer_type, y = Mean, label = paste('Mean:',Mean)), size = 4, color='black') + 
  geom_text(data = stats, aes(x = primer_type, y = Max*1.01, label = paste('Max:',Max)), size = 4, color='black') + 
  labs(x='Primer Type', y='Primer Size') + 
  scale_color_brewer(palette="Set2") + scale_fill_brewer(palette="Set2") +
  theme(axis.text=element_text(size=10), axis.title=element_text(size=10), 
        plot.title=element_text(size=15), plot.caption=element_text(size=10), 
        legend.title=element_text(size=10), legend.text=element_text(size=10))
dev.off()



# GC dist

mins <- df %>% group_by(primer_type) %>% summarise_at(vars(pct_gc), list(name = min))
means <- df %>% group_by(primer_type) %>% summarise_at(vars(pct_gc), list(name = mean))
medians <- df %>% group_by(primer_type) %>% summarise_at(vars(pct_gc), list(name = median))
maxs <- df %>% group_by(primer_type) %>% summarise_at(vars(pct_gc), list(name = max))
stats <- data.frame(mins, medians, means, maxs)[,c(1,2,4,6,8)]
colnames(stats) <- c('primer_type', 'Min', 'Median','Mean', 'Max')
stats <- stats %>% mutate_if(is.numeric, ~ . * 100)

pdf(paste0(out,'.primer_gc.pdf'), width = 6, height = 6)
ggplot(df, aes(x=primer_type, y=pct_gc*100, color=primer_type, fill=primer_type)) +
  geom_boxplot(alpha=0.5, outlier.color=NA) + geom_jitter(shape=16, width=0.2, height=0, size=2) + 
  geom_text(data = stats, aes(x = primer_type, y = Min, label = paste('Min:',Min)), size = 4, color='black') + 
  geom_text(data = stats, aes(x = primer_type, y = Mean, label = paste('Mean:',Mean)), size = 4, color='black') + 
  geom_text(data = stats, aes(x = primer_type, y = Max*1.01, label = paste('Max:',Max)), size = 4, color='black') + 
  labs(x='Primer Type', y='Primer GC%') + 
  scale_color_brewer(palette="Set2") + scale_fill_brewer(palette="Set2") +
  theme(axis.text=element_text(size=10), axis.title=element_text(size=10), 
        plot.title=element_text(size=15), plot.caption=element_text(size=10), 
        legend.title=element_text(size=10), legend.text=element_text(size=10))
dev.off()



