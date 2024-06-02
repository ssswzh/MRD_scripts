# path
# /public2/home/lizhe/zhangsw/jinyun/project/PRJZ_001_Pan_gynecologic_wes_mutation_research/Standard_Analysis/Chemotherapy/
# /public2/home/lizhe/zhangsw/jinyun/project/PRJZ_001_Pan_gynecologic_wes_mutation_research/Standard_Analysis/PAPRI/20210830_WES

# jupiter
# /mnt/ddngs/zhangshouwei/project/PRJZ_001_Pan_gynecologic_wes_mutation_research/PAPRI/Basis_results/Results/MC3

# <<<--- variables --->>>

tumorBam='P13.tumor.bam'
tumorSampleName='P13_T'
normalBam='P13.normal.bam'
normalSampleName='P13_N'
out='P13'
outdir='/mnt/ddngs/zhangsw/project/MRD/OV_WES/P13'
#gatk4='/mnt/ddngs/zhangsw/bin/gatk-package-4.1.0.0-local.jar'
gatk4='/mnt/ddngs/zhangsw/miniconda3/bin/gatk'
refFa='/mnt/user/errand/pipeline/MC3/ref/Homo_sapiens_assembly19.fasta'
bed='/mnt/ddngs/zhangsw/project/MRD/OV_WES/bed/Agilent_v6r2.Covered.bed'
scriptDir='/mnt/ddngs/zhangsw/project/MRD/OV_WES/scripts'
annovarDB='/resource/annovar/humandb'
annovarscr='/mnt/ddngs/zhangsw/software/annovar'

# <<<--- SNV indel calling --->>>

## mutect

mkdir -p ${outdir}/mutect2
${gatk4} --java-options "-XX:ParallelGCThreads=4 -Xmx16G " Mutect2 \
    -R ${refFa} \
    -I ${tumorBam} -I ${normalBam} \
    -normal ${normalSampleName} -tumor ${tumorSampleName} \
    -L ${bed} \
    --native-pair-hmm-threads 8 \
    --f1r2-tar-gz ${outdir}/mutect2/${out}.mutect2F1R2.tar.gz \
    --max-reads-per-alignment-start 0 \
    --initial-tumor-lod 0 \
    -O ${outdir}/mutect2/${out}.mutect2raw.vcf
    #--germline-resource /resource/annovar/humandb/hg19_gnomad211_exome.txt 

${gatk4} --java-options "-XX:ParallelGCThreads=4 -Xmx8G " LearnReadOrientationModel \
    -I ${outdir}/mutect2/${out}.mutect2F1R2.tar.gz \
    -O ${outdir}/mutect2/${out}.mutect2ReadOrientation.tar.gz

${gatk4} --java-options "-XX:ParallelGCThreads=4 -Xmx8G " GetPileupSummaries \
    -I ${tumorBam} \
    -V /mnt/ddngs/zhangsw/database/GRCh37/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz \
    -L ${bed} \
    -O ${outdir}/mutect2/${out}.mutect2GetPileupSummaries.table

${gatk4} --java-options "-XX:ParallelGCThreads=4 -Xmx8G " CalculateContamination \
        -I ${outdir}/mutect2/${out}.mutect2GetPileupSummaries.table \
        -tumor-segmentation ${outdir}/mutect2/${out}.mutect2Segments.table \
        -O ${outdir}/mutect2/${out}.mutect2CalculateContamination.table

${gatk4} --java-options "-XX:ParallelGCThreads=4 -Xmx8G " FilterMutectCalls \
    -R ${refFa} \
    -V ${outdir}/mutect2/${out}.mutect2raw.vcf \
    -O ${outdir}/mutect2/${out}.mutect2Filter.vcf \
    --stats ${outdir}/mutect2/${out}.mutect2raw.vcf.stats \
    --max-events-in-region 10 \
    --min-median-mapping-quality 20 \
    --tumor-segmentation ${outdir}/mutect2/${out}.mutect2Segments.table \
    --contamination-table ${outdir}/mutect2/${out}.mutect2CalculateContamination.table \
    --ob-priors ${outdir}/mutect2/${out}.mutect2ReadOrientation.tar.gz

### filter AD DP AF, different for snv and indel

bcftools norm ${outdir}/mutect2/${out}.mutect2Filter.vcf -m - -Ov -o ${outdir}/mutect2/${out}.mutect2Filter.norm.vcf
bgzip --force ${outdir}/mutect2/${out}.mutect2Filter.norm.vcf && tabix ${outdir}/mutect2/${out}.mutect2Filter.norm.vcf.gz
bcftools view -i 'FILTER=="PASS" && (FORMAT/AD[1:1]/FORMAT/DP[1]) >= 0.01 && FORMAT/DP[1] >= 30 && FORMAT/AD[0:1] <= 5 && (FORMAT/AD[0:1]/FORMAT/DP[0]) <= 0.01 && FORMAT/DP[0] >= 10' -Ov -o ${outdir}/mutect2/${out}.mutect2Filter.norm.snp.vcf ${outdir}/mutect2/${out}.mutect2Filter.norm.vcf.gz
bcftools view -i 'FILTER=="PASS" && (FORMAT/AD[1:1]/FORMAT/DP[1]) < 0.02 || FORMAT/DP[1] <=50' -Ov ${outdir}/mutect2/${out}.mutect2Filter.norm.vcf.gz | grep -v '#' |awk -F '\t' 'length($4)!=length($5)' > ${outdir}/mutect2/${out}.mutect2Filter.norm.indelFailed.vcf
grep -v -f ${outdir}/mutect2/${out}.mutect2Filter.norm.indelFailed.vcf ${outdir}/mutect2/${out}.mutect2Filter.norm.snp.vcf > ${outdir}/mutect2/${out}.mutect2Filter.norm.snp.filtered.vcf


## varscan

### pre-process

mkdir -p ${outdir}/varscan2somatic
samtools mpileup \
    --min-MQ 20 --min-BQ 20 --adjust-MQ 50 --no-BAQ \
    --positions ${bed} --fasta-ref ${refFa} \
    --output ${out}.normal.pileup ${normalBam}

samtools mpileup \
    --min-MQ 20 --min-BQ 20 --adjust-MQ 50 --no-BAQ \
    --positions ${bed} --fasta-ref ${refFa} \
    --output ${out}.tumor.pileup ${tumorBam}

### varscan somatic calling. The command will generate two output files, one for SNVs (output.basename.snp) and one for indels (output.basename.indel). 

varscan somatic \
    ${out}.normal.pileup ${out}.tumor.pileup ${outdir}/varscan2somatic/${out}.varscanSomatic \
    --min-coverage-normal 10 --min-var-freq 0.01 --tumor-purity 0.5 --output-vcf 1

varscan processSomatic \
    ${outdir}/varscan2somatic/${out}.varscanSomatic.snp.vcf \
    --min-tumor-freq 0.01 --max-normal-freq 0.01 --p-value 0.01
varscan processSomatic \
    ${outdir}/varscan2somatic/${out}.varscanSomatic.indel.vcf \
    --min-tumor-freq 0.01 --max-normal-freq 0.01 --p-value 0.001

### bam-readcount filter SNV, change min-var-frac=0.01 in fpfilter.pl

#cat ${outdir}/varscan2somatic/${out}.varscanSomatic.snp.Somatic.vcf| vcf2bed - > ${outdir}/varscan2somatic/${out}.varscanSomatic.snp.Somatic.bed
/mnt/ddngs/zhangsw/miniconda3/envs/bam-readcount/bin/bam-readcount \
    --min-mapping-quality 1 --min-base-quality 20 --max-warnings 0 \
    -f ${refFa} -l ${bed} \
    ${tumorBam} > ${outdir}/varscan2somatic/${out}.varscanSomatic.snp.Somatic.bam-readcount

perl ${scriptDir}/fpfilter.pl \
    --snp-file ${outdir}/varscan2somatic/${out}.varscanSomatic.snp.Somatic.vcf --readcount-file ${outdir}/varscan2somatic/${out}.varscanSomatic.snp.Somatic.bam-readcount

### varscan filter SNV and indel, only keep high-confidence indels

varscan somaticFilter ${outdir}/varscan2somatic/${out}.varscanSomatic.snp.Somatic.fp_pass \
    --min-coverage 30 --min-reads2 5 --min-var-freq 0.01 --p-value 0.01 \
    --indel-file ${outdir}/varscan2somatic/${out}.varscanSomatic.indel \
    --output-file ${outdir}/varscan2somatic/${out}.varscanSomatic.snp.Somatic.fp_pass.somaticFilter.vcf
varscan somaticFilter ${outdir}/varscan2somatic/${out}.varscanSomatic.indel.Somatic.hc.vcf \
    --min-coverage 50 --min-reads2 10 --min-var-freq 0.01 --p-value 0.001 \
    --output-file ${outdir}/varscan2somatic/${out}.varscanSomatic.indel.Somatic.hc.somaticFilter.vcf

### varscan merge SNV and indels

bgzip --force ${outdir}/varscan2somatic/${out}.varscanSomatic.snp.Somatic.fp_pass.somaticFilter.vcf && tabix ${outdir}/varscan2somatic/${out}.varscanSomatic.snp.Somatic.fp_pass.somaticFilter.vcf.gz
bgzip --force ${outdir}/varscan2somatic/${out}.varscanSomatic.indel.Somatic.hc.somaticFilter.vcf && tabix ${outdir}/varscan2somatic/${out}.varscanSomatic.indel.Somatic.hc.somaticFilter.vcf.gz

bcftools concat \
    -a -R ${bed} -O v ${outdir}/varscan2somatic/${out}.varscanSomatic.snp.Somatic.fp_pass.somaticFilter.vcf.gz ${outdir}/varscan2somatic/${out}.varscanSomatic.indel.Somatic.hc.somaticFilter.vcf.gz|\
    bcftools sort --temp-dir ./ - > ${outdir}/varscan2somatic/${out}.varscanSomatic.totalmutation.vcf


## merge varscan and mutect2, filtering, get true positives

### merge, filter AD DP

mkdir -p ${outdir}/mergeMutations
bgzip --force ${outdir}/varscan2somatic/${out}.varscanSomatic.totalmutation.vcf && tabix ${outdir}/varscan2somatic/${out}.varscanSomatic.totalmutation.vcf.gz
bgzip --force ${outdir}/mutect2/${out}.mutect2Filter.norm.snp.filtered.vcf && tabix ${outdir}/mutect2/${out}.mutect2Filter.norm.snp.filtered.vcf.gz

bcftools merge \
    --force-samples -R ${bed} -O v ${outdir}/varscan2somatic/${out}.varscanSomatic.totalmutation.vcf.gz ${outdir}/mutect2/${out}.mutect2Filter.norm.snp.filtered.vcf.gz| \
    bcftools sort --temp-dir ./ - > ${outdir}/mergeMutations/${out}.merged.vcf

bgzip --force ${outdir}/mergeMutations/${out}.merged.vcf && tabix ${outdir}/mergeMutations/${out}.merged.vcf.gz
bcftools view -i '(FORMAT/FREQ[1:0] >= 0.01 && FORMAT/AF[3:0] >= 0.01) || (FORMAT/FREQ[1:0] >= 0.05)' -Ov -o ${outdir}/mergeMutations/${out}.merged.filterAF.vcf ${outdir}/mergeMutations/${out}.merged.vcf.gz

### filter blacklist (Duke, DAC, repeatMasker, simpleRepeat)

bedtools subtract -nonamecheck -header -a ${outdir}/mergeMutations/${out}.merged.filterAF.vcf -b /mnt/ddngs/zhangsw/database/GRCh37/wgEncodeDacMapabilityConsensusExcludable.bed.gz |\
    bedtools subtract -nonamecheck -header -a - -b /mnt/ddngs/zhangsw/database/GRCh37/wgEncodeDukeMapabilityRegionsExcludable.bed.gz |\
    bedtools subtract -nonamecheck -header -a - -b /mnt/ddngs/zhangsw/database/GRCh37/hg19_repeatMasker.bed |\
    bedtools subtract -nonamecheck -header -a - -b /mnt/ddngs/zhangsw/database/GRCh37/simpleRepeat.bed > ${outdir}/mergeMutations/${out}.merged.filterAF.filterBlacklist.vcf


## annotation

### annovar

perl ${annovarscr}/table_annovar.pl ${outdir}/mergeMutations/${out}.merged.filterAF.filterBlacklist.vcf ${annovarDB} --remove \
    --vcfinput -nastring . --otherinfo  --buildver hg19 \
    --protocol  refGene,ljb26_all,clinvar_20200316,intervar_20180118,cosmic90,EAS.sites.2015_08,esp6500siv2_all,exac03,avsnp150,gnomad211_exome \
    --operation g,f,f,f,f,f,f,f,f,f \
    --outfile ${outdir}/mergeMutations/${out}.annovar

### vep annotation

#docker run --rm -u 1925:1000 -v /mnt/ddngs/zhangsw/database/vep104/:/database/ -v /mnt/user/errand/pipeline/MC3/ref/:/ref/ -v ${outdir}/mergeMutations/:/outdir/ registry.cn-beijing.aliyuncs.com/base-docker/ensembl-vep:release_104.3 /opt/vep/src/ensembl-vep/vep cache \
    --offline --no_progress -no_stats \
    --hgvs --pubmed --force_overwrite --exclude_predicted \
    --symbol --refseq --format vcf  --distance 0 \
    --dir_cache /database/ --fasta /ref/Homo_sapiens_assembly19.fasta \
    --input_file /outdir/${out}.merged.filterAF.filterBlacklist.vcf --output_file /outdir/${out}.vep.vcf \
    --vcf --af_1kg --af_gnomad --no_escape

### merge annotations

#python ${scriptDir}/annotation_combine.weib.py \
    --annovar ${outdir}/mergeMutations/${out}.annovar.hg19_multianno.txt \
    --vep ${outdir}/mergeMutations/${out}.vep.vcf \
    --vcf_sample_name Normal,Tumor,${out}_N,${out}_T \
    --select_transcript_file /mnt/ddngs/zhangsw/database/GRCh37/GCF_000001405.25_GRCh37.p13_genomic.selected_transcript.txt \
    --combine_strategy union \
    --output ${outdir}/mergeMutations/${out}.annotationCombined.txt




# <<<--- Copy Number Variation --->>>

## VarScan copynumber

mkdir -p ${outdir}/varscan2copynumber
varscan copynumber \
    ${out}.normal.pileup ${out}.tumor.pileup ${outdir}/varscan2copynumber/${out}.varscanCopynumber \
    --min-coverage 8 --min-segment-size 50 --data-ratio 0.6 
varscan copyCaller \
    ${outdir}/varscan2copynumber/${out}.varscanCopynumber.copynumber \
    --output-file ${outdir}/varscan2copynumber/${out}.varscanCopynumber.called \
    --homdel-file ${outdir}/varscan2copynumber/${out}.varscanCopynumber.homdel \
    --min-coverage 8 


## platypus get germline variants, calculate BAF

/mnt/ddngs/zhangsw/miniconda3/envs/platypus/bin/platypus callVariants \
    --genIndels 0 \
    --refFile ${refFa} \
    --bamFiles ${normalBam} \
    --regions ${bed} \
    --minReads 20 \
    --output ${outdir}/varscan2copynumber/${out}.normal.platypus.vcf

python ${scriptDir}/calculate_baf_platypus.py \
    --vcf ${outdir}/varscan2copynumber/${out}.normal.platypus.vcf \
    --out ${outdir}/varscan2copynumber/${out}.normal.platypus.baf \
    --fmt both

## extract somatic BAF from germline tsv

python ${scriptDir}/somatic_baf_by_bamcount.py \
    --germline ${outdir}/varscan2copynumber/${out}.normal.platypus.baf.tsv \
    --somatic ${outdir}/varscan2somatic/${out}.varscanSomatic.snp.Somatic.bam-readcount \
    --out ${outdir}/varscan2copynumber/${out}.tumor.platypus.baf.tsv

## extract logR from BAF file

python ${scriptDir}/extract_mutation_cnv.py \
    --baf ${outdir}/varscan2copynumber/${out}.normal.platypus.baf.tsv \
    --cnv ${outdir}/varscan2copynumber/${out}.varscanCopynumber.called \
    --out ${outdir}/varscan2copynumber/${out}.varscanCopynumber.called.logR.tsv
awk 'BEGIN{FS="\t";OFS="\t"}NR==1{print $0}NR>1{print $1,$2,$3,0}' ${outdir}/varscan2copynumber/${out}.varscanCopynumber.called.logR.tsv > ${outdir}/varscan2copynumber/${out}.normal.logR.tsv

## run ASCAT

Rscript --vanilla ${scriptDir}/ascat_cnv.r \
    --somaticcnv ${outdir}/varscan2copynumber/${out}.varscanCopynumber.called.logR.tsv \
    --somaticbaf ${outdir}/varscan2copynumber/${out}.tumor.platypus.baf.tsv \
    --germlinecnv ${outdir}/varscan2copynumber/${out}.normal.logR.tsv \
    --germlinebaf ${outdir}/varscan2copynumber/${out}.normal.platypus.baf.tsv \
    --outdir ${outdir}/varscan2copynumber/ \
    --out ${out}.varscanCopynumber.called.ascat



# <<<--- PyClone cluster and phylogenetic tree construction --->>>

## prepare data

mkdir -p ${outdir}/pyclone_analysis/
python ${scriptDir}/pyclone_dat_generation_vcf.py \
    --vcf ${outdir}/mergeMutations/${out}.merged.filterAF.filterBlacklist.vcf \
    --cnv ${outdir}/varscan2copynumber/${out}.varscanCopynumber.called.ascat \
    --out ${outdir}/pyclone_analysis/${out}.pyclone_dat.tsv \
    --gender female --sampleid ${out} --depth 30 --ad 5 --af 0.01

## run new pyclone

/mnt/ddngs/zhangsw/miniconda3/envs/pyclone-vi/bin/pyclone-vi fit \
    -i ${outdir}/pyclone_analysis/${out}.pyclone_dat.tsv \
    -o ${outdir}/pyclone_analysis/${out}.pyclone_dat.h5 \
    --num-clusters 40 --num-restarts 10 --max-iters 10000 \
    --density binomial --seed 93521 > ${outdir}/pyclone_analysis/log.pyclone 2>&1
/mnt/ddngs/zhangsw/miniconda3/envs/pyclone-vi/bin/pyclone-vi write-results-file \
    -i ${outdir}/pyclone_analysis/${out}.pyclone_dat.h5 \
    -o ${outdir}/pyclone_analysis/${out}.pyclone_dat.h5.tsv

# old pyclone
#sed -i 's/alt_counts/var_counts/' ${outdir}/pyclone_analysis/${out}.pyclone_dat.tsv
#/mnt/ddngs/zhangsw/miniconda3/envs/pyclone/bin/PyClone run_analysis_pipeline \
    --in_files ${outdir}/pyclone_analysis/${out}.pyclone_dat.tsv \
    --working_dir ${outdir}/pyclone_analysis/ \
    --burnin 1000 --num_iters 10000 --max_clusters 40 \
    --density pyclone_binomial --seed 93521 > ${outdir}/pyclone_analysis/log.pyclone 2>&1

## citup phylogenetic tree construction

#sort -k1 P13.pyclone_dat.h5.tsv|cut -f3|grep -v "cluster_id" > test.citup_cluster.txt
#sort -k1 P13.pyclone_dat.tsv|cut -f5|grep -v "vaf" > test.citup_freq.txt
#/mnt/ddngs/zhangsw/miniconda3/envs/citup/bin/run_citup_qip.py --submit local --tmpdir citup_tmp.pyclone-vi test.citup_freq.txt test.pyclone-vi.citup_cluster.txt test.pyclone-vi.citup_cluster.h5
