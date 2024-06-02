library(pctGCdata, lib="/mnt/user/errand/pipeline/OncoDeficiencyPro/tools/R")
library(facets, lib="/mnt/user/errand/pipeline/OncoDeficiencyPro/tools/R")


args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3)
    stop("Usage: Rscript run_facet.R <snp-fileup_file> <sample_name> <outdir>")

infile <- args[1]
sample <- args[2]
outdir <- args[3]
dir.create(outdir)
set.seed(1234)

rcmat <- readSnpMatrix(infile)
write.table(rcmat, file=paste(outdir, "/", sample, ".1.site.txt", sep = ''), quote = F, sep = "\t", row.names=F)

xx <- preProcSample(rcmat)
write.table(xx$jointseg, file=paste(outdir, "/", sample, ".2.site.segment.txt", sep = ''), quote = F, sep = "\t", row.names = F)
#如果segment过多，需要调参数
oo <- procSample(xx, cval = 150)
write.table(oo$out, file=paste(outdir, "/", sample, ".3.segment.cluster.cf.cn.txt", sep = ''), quote = F, sep = "\t", row.names = F)
write.table(oo$jointseg, file=paste(outdir, "/", sample, ".3.site.segment.cluster.txt", sep = ''), quote = F, sep = "\t", row.names = F)
write.table(oo$dipLogR, file=paste(outdir, "/", sample, ".3.dipLogR.txt", sep = ''), quote = F, sep = "\t", row.names = F)
write.table(oo$alBalLogR, file=paste(outdir, "/", sample, ".3.alBalLogR.txt", sep = ''), quote = F, sep = "\t", row.names = F)
write.table(oo$mafR.thresh, file=paste(outdir, "/", sample, ".3.mafR.thresh.txt", sep = ''), quote = F, sep = "\t", row.names = F)
write.table(oo$flags, file=paste(outdir, "/", sample, ".3.flags1.txt", sep = ''), quote = F, sep = "\t", row.names = F)



fit <- emcncf(oo)
# zsw
fit$cncf$mcn.em <- fit$cncf$tcn.em - fit$cncf$lcn.em
if(is.na(fit$purity)) {
  fit$cncf$purity <- 1
}else{
  fit$cncf$purity <- fit$purity
}
fit$cncf$ploidy <- fit$ploidy
write.table(fit$cncf, file=paste(outdir, "/", sample, ".4.segment.final.txt", sep = ''), quote = F, sep = "\t", row.names = F)
#table <- data.frame(sampleid=sample, loglik=fit$loglik, purity=fit$purity, ploidy=fit$ploidy, dipLogR=fit$dipLogR)
#write.table(table, file=paste(outdir, "/", sample, ".4.pp.final.txt", sep = ''), quote = F, sep = "\t", row.names = F)


save(oo, fit, file = file.path(outdir, paste(sample, '.rda', sep = '')))

pdf(file.path(outdir, paste(sample, "_seg.pdf", sep = '')),
    width = 9, height = 7)
plotSample(x = oo, emfit = fit, plot.type = 'both', sname = sample)
dev.off()

pdf(file.path(outdir, paste(sample, "_diag.pdf", sep = '')),
    width = 7, height = 7)
logRlogORspider(oo$out, oo$dipLogR)
dev.off()
