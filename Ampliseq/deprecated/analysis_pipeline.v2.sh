#!/usr/bin/bash
# -*- coding: utf-8 -*-
# @Time    : 2022/08/01
# @Author  : zhangsiwen
# @contact : zhang.siwen
# @LastModified: 2022/08/01
# @ChangeLog
#     20220801, first version
#     20220919, add Rscript boxplot_with_quantiles.r to draw target coverage boxplot
#     20221008, third version, add fastp


# <<<--- variables --->>>

shellname=$0
scriptDir=`dirname ${shellname}`
VardictBin='/mnt/ddngs/zhangsw/software/VarDict-1.8.3/bin'
#refFa='/mnt/user/errand/pipeline/MC3/ref/Homo_sapiens_assembly19.fasta'
refFa='/mnt/ddngs/zhangsw/database/GRCh37_clinical/nochrY/Homo_sapiens_assembly19.nochrY.fasta'
Rscript='/mnt/ddngs/zhangsw/miniconda3/bin/Rscript'
outdir=`pwd`'/outdir'
prefix='sample'
thread=4

# usage
function usage {
    printf -- "Usage:\n    sh ${shellname} -r1 R1.fq -r2 R2.fq -l bed -s site.bed -v refFa -o outdir -p prefix -t thread \n\n";
    printf -- "Required arguments: \n";
    printf -- "    -f|--forward       R1.fq \n";
    printf -- "    -r|--reverse       R2.fq \n";
    printf -- "    -l|--bed           bed region \n";
    printf -- "    -s|--site          site coordinate bed \n";
    printf -- "Optional arguments: \n";
    printf -- "    -v|--refFa         reference fasta, default '${refFa}' \n"
    printf -- "    -o|--outdir        output path, default '${outdir}' \n";
    printf -- "    -p|--prefix        output prefix, default '${prefix}' \n";
    printf -- "    -t|--thread        thread, default '${thread}' \n";
    printf -- "    -h|--help          help \n";
    printf -- "@Author: Siwen Zhang \n\n";
    exit 0
}


# Arguments passing
TEMP=`getopt -o f:r:l:s:v:o:p:t:h --long forward:,reverse:,bed:,site:,refFa:,outdir:,prefix:,thread:,help \
     -n ${shellname} -- "$@"`

# no arguments detected
if [ $? != 0 ] || [ $# == 0 ]
    then usage
fi

# set arguments
eval set -- "$TEMP"
while true ; do
    case "$1" in
        -f|--forward)
            forward="$2"; shift 2 ;;
        -r|--reverse)
            reverse="$2"; shift 2 ;;
        -l|--bed)
            bed="$2"; shift 2 ;;
        -s|--site)
            site="$2"; shift 2 ;;
        -v|--refFa)
            refFa="$2"; shift 2 ;;
        -o|--outdir)
            outdir="$2"; shift 2 ;;
        -p|--prefix)
            prefix="$2"; shift 2 ;;
        -t|--thread)
            thread="$2"; shift 2 ;;
        -h|--help)
            usage ;;
        --)
            shift ; break; exit 1;;
        ?)
            usage ;;
        *)
            usage ;;
    esac
done

# arguments not null
if [ -f ${forward} ] && [ -f ${reverse} ] && [ -f ${bed} ] && [ -f ${site} ]
then echo "Files found."
    echo "R1: ${forward}"
    echo "R2: ${reverse}"
    echo "BED: ${bed}"
    echo "SITE: ${site}"
else usage
fi








# pear merge fq1 and fq2

mkdir -p ${outdir}/fastq_merge 

# usearch merge fq1 and fq2
#https://drive5.com/usearch/manual/cmds_all.html
#https://drive5.com/usearch/manual/exp_errs.html

# debug mode add:  -alnout ${outdir}/fastq_merge/${prefix}.merged.aln 

if [[ ${forward} == *.gz ]]
then 
    gunzip -c ${forward} > ${outdir}/fastq_merge/R1.fq; forward=${outdir}/fastq_merge/R1.fq
    gunzip -c ${reverse} > ${outdir}/fastq_merge/R2.fq; reverse=${outdir}/fastq_merge/R2.fq
fi

usearch -fastq_mergepairs ${forward} \
    -reverse ${reverse} \
    -fastqout ${outdir}/fastq_merge/${prefix}.merged.fq \
    -tabbedout ${outdir}/fastq_merge/${prefix}.merged.tsv \
    -alnout ${outdir}/fastq_merge/${prefix}.merged.aln \
    -fastqout_notmerged_fwd ${outdir}/fastq_merge/${prefix}.notmerged_R1.fq \
    -fastqout_notmerged_rev ${outdir}/fastq_merge/${prefix}.notmerged_R2.fq \
    -threads ${thread} \
    -fastq_maxdiffs 10 -fastq_pctid 90 -fastq_merge_maxee 0.01 \
    -fastq_minovlen 40 -fastq_minmergelen 40 -fastq_maxmergelen 200 \
    2> ${outdir}/fastq_merge/${prefix}.merged.stats

sed -n '10,11p' ${outdir}/fastq_merge/${prefix}.merged.stats|awk 'NR==1{total=$1;print "Pair_Reads\t"$1}NR==2{n=$1;print "Merged_Reads\t"$1}END{print "Merged_Reads_Pct\t"n/total}' > ${outdir}/fastq_merge/${prefix}.merged.stats.tsv

if [ -f ${outdir}/fastq_merge/R1.fq ]
then rm ${outdir}/fastq_merge/R1.fq ${outdir}/fastq_merge/R2.fq
fi


# ----- merged reads ----- #

# mapping

mkdir -p ${outdir}/mapping

bwa mem -t ${thread} -M \
    -R "@RG\\tID:NULL\\tSM:${prefix}\\tLB:NULL\\tPU:NULL\\tPL:illumina\\tCN:prec_sci" \
    ${refFa} ${outdir}/fastq_merge/${prefix}.merged.fq 2> log.txt 1> ${outdir}/mapping/${prefix}.sam

samtools view -Sb -@ ${thread} -t ${refFa}.fai ${outdir}/mapping/${prefix}.sam > ${outdir}/mapping/${prefix}.bam

sambamba sort -t ${thread} -l 1 -o ${outdir}/mapping/${prefix}.sorted.bam ${outdir}/mapping/${prefix}.bam

samtools index -@ ${thread} ${outdir}/mapping/${prefix}.sorted.bam ${outdir}/mapping/${prefix}.sorted.bam.bai 


# stats

mkdir -p ${outdir}/stats

# count base number for each region
#https://broadinstitute.github.io/picard/

if [ -f ${refFa}.dict ]
then refDict=${refFa}.dict
else refDict=`echo ${refFa}|sed 's/fasta$//'|sed 's/fa$//'`"dict"
fi

bed_name=`basename -s '.bed' ${bed}`

picard BedToIntervalList \
    -I ${bed} \
    -O ${outdir}/stats/${bed_name}.intervalList \
    -SD ${refDict}

picard CollectTargetedPcrMetrics \
    I=${outdir}/mapping/${prefix}.sorted.bam \
    O=${outdir}/stats/${prefix}.CollectTargetedPcrMetrics \
    R=${refFa} \
    AMPLICON_INTERVALS=${outdir}/stats/${bed_name}.intervalList \
    TARGET_INTERVALS=${outdir}/stats/${bed_name}.intervalList \
    PER_TARGET_COVERAGE=${outdir}/stats/${prefix}.per_target_cov \
    PER_BASE_COVERAGE=${outdir}/stats/${prefix}.per_base_cov \
    MINIMUM_BASE_QUALITY=20 \
    NEAR_DISTANCE=100 \
    COVERAGE_CAP=1000000 
    #MINIMUM_MAPPING_QUALITY=30 

python ${scriptDir}/mapping_qc_stats.py \
    --metrics ${outdir}/stats/${prefix}.CollectTargetedPcrMetrics \
    --per_base_cov ${outdir}/stats/${prefix}.per_base_cov \
    --per_target_cov ${outdir}/stats/${prefix}.per_target_cov \
    --bed ${bed} \
    --out ${outdir}/stats/${prefix}.QCstats.tsv

${Rscript} --vanilla ${scriptDir}/boxplot_with_quantiles.r \
    --in ${outdir}/stats/${prefix}.per_target_cov \
    --out ${prefix}.per_target_cov \
    --dir ${outdir}/stats/ \
    --key mean_coverage


# mutations

mkdir ${outdir}/mutations

# pileup-bcftools

#samtools mpileup -t DP -t SP -t AD --VCF --uncompressed --max-depth 500000 \
#    -o ${outdir}/mutations/${prefix}.mpileup.vcf \
#    -f ${refFa} -l ${bed} --min-MQ 30 --min-BQ 20 \
#    ${outdir}/mapping/${prefix}.sorted.bam
#bcftools call -v -m -P 0.0001 -T ${bed} -A ${outdir}/mutations/${prefix}.mpileup.vcf|bcftools norm -O v -m - - > ${outdir}/mutations/${prefix}.mpileup.filter.vcf

# vardict

awk 'BEGIN{FS="\t";OFS="\t"}{if($1!=""){print $1,$2-20,$3+20,$4,"0",".",$2,$3}}' ${bed} > ${outdir}/mutations/${bed_name}.amplicon.bed

${VardictBin}/VarDict -G ${refFa} \
    -N "SAMPLE" \
    -b "${outdir}/mapping/${prefix}.sorted.bam" \
    -th ${thread} -O 10 -Q 30 -F 0 -p -a 20:0.2 \
    --nosv -f 0.00001 -c 1 -S 2 -E 3 -g 4 ${outdir}/mutations/${bed_name}.amplicon.bed | \
    ${VardictBin}/teststrandbias.R | \
    ${VardictBin}/var2vcf_valid.pl \
        -A -N "SAMPLE" -d 1000 -v 1 -f 0.00001 -P 0 > ${outdir}/mutations/${prefix}.vardict.vcf

# amplicon mode bed file. There first two numbers include primers, and the last two numbers are insert only.
# VarDict will perform in silico trimming of primers, and for variants in overlapping amplicons, they needs to be detectable by both amplicons, otherwise, they're flagged as amplicon-bias, and likely false positives.
#chr1    115247094       115247253       NRAS    0       .       115247117       115247232
#chr1    115247202       115247341       NRAS    0       .       115247224       115247323

# filter

bgzip -f ${outdir}/mutations/${prefix}.vardict.vcf && tabix ${outdir}/mutations/${prefix}.vardict.vcf.gz
bcftools view -i "FORMAT/DP>1000 && FORMAT/VD>=10 && %FILTER='PASS'" ${outdir}/mutations/${prefix}.vardict.vcf.gz > ${outdir}/mutations/${prefix}.vardict.filter.vcf


awk 'BEGIN{FS="\t";OFS="\t"}
    NR==FNR{if($1!~"#"){split($10,info,":");af=info[3]/info[2]; a[$1":"$2":"$4":"$5]=$4">"$5"\t"info[2]"\t"info[3]"\t"af}}
    NR>FNR{if($1=="#Chromosome"){print $0,"detected","depth","reads","vaf"} else{if(a[$1":"$3":"$4":"$5]){print $0,a[$1":"$3":"$4":"$5]}else{print $0,$4">"$5,"0","0","0"}}}' ${outdir}/mutations/${prefix}.vardict.filter.vcf ${site}| \
    awk 'BEGIN{FS="\t";OFS="\t"}
    NR==FNR{a[$1":"$2]=$4}
    NR>FNR{if($10==0){print $1,$2,$3,$4,$5,$6,$7,$8,$4">"$5,a[$1":"$3],"0","0"} else{print $0}}' ${outdir}/stats/${prefix}.per_base_cov - > ${outdir}/mutations/${prefix}.vardict.filter.intersect.tsv

awk 'BEGIN{FS="\t";OFS="\t"}NR==1{print $0,"result"}NR>1{if($NF>=0.0003){print $0,"Positive"}else{print $0,"Negative"}}' ${outdir}/mutations/${prefix}.vardict.filter.intersect.tsv > ${outdir}/mutations/${prefix}.vardict.filter.intersect.result.tsv


#awk 'BEGIN{FS="\t";OFS="\t"} NR==FNR{split($10,info,":"); a[$1":"$2":"$4":"$5]=info[3]/info[2]} NR>FNR{print $0,a[$1":"$3":"$4":"$5]}' ${outdir}/mutations/${prefix}.vardict.filter.intersect.tsv ${site} > ${outdir}/mutations/${prefix}.vardict.filter.intersect.only_vaf.tsv

#awk '{if($1~/^#/){print $0} else{split($10,a,":"); if(($4=="G"&&$5=="A"&&a[3]/a[2]>=0.0003)||($4=="C"&&$5=="T"&&a[3]/a[2]>=0.0003)){print $0}}}' ${outdir}/mutations/${prefix}.vardict.filter.intersect.tsv > ${outdir}/mutations/${prefix}.vardict.filter.intersect.filter_transition.tsv





# ----- not merged reads ----- #

mkdir -p ${outdir}/not_merged

# mapping

bwa mem -t ${thread} -M -R "@RG\\tID:NULL\\tSM:${prefix}\\tLB:NULL\\tPU:NULL\\tPL:illumina\\tCN:prec_sci" ${refFa} ${outdir}/fastq_merge/${prefix}.notmerged_R1.fq ${outdir}/fastq_merge/${prefix}.notmerged_R2.fq | \
    samtools view -Sb -@ ${thread} -t ${refFa}.fai - > ${outdir}/not_merged/${prefix}.bam

sambamba sort -t ${thread} -l 1 -o ${outdir}/not_merged/${prefix}.sorted.bam ${outdir}/not_merged/${prefix}.bam

samtools index -@ ${thread} ${outdir}/not_merged/${prefix}.sorted.bam ${outdir}/not_merged/${prefix}.sorted.bam.bai


# stats

picard CollectTargetedPcrMetrics \
    I=${outdir}/not_merged/${prefix}.bam \
    O=${outdir}/not_merged/${prefix}.CollectTargetedPcrMetrics \
    R=${refFa} \
    AMPLICON_INTERVALS=${outdir}/stats/${bed_name}.intervalList \
    TARGET_INTERVALS=${outdir}/stats/${bed_name}.intervalList \
    PER_TARGET_COVERAGE=${outdir}/not_merged/${prefix}.per_target_cov \
    PER_BASE_COVERAGE=${outdir}/not_merged/${prefix}.per_base_cov \
    MINIMUM_BASE_QUALITY=20 \
    MINIMUM_MAPPING_QUALITY=30 \
    NEAR_DISTANCE=100 

python ${scriptDir}/mapping_qc_stats.py \
    --metrics ${outdir}/not_merged/${prefix}.CollectTargetedPcrMetrics \
    --per_base_cov ${outdir}/not_merged/${prefix}.per_base_cov \
    --per_target_cov ${outdir}/not_merged/${prefix}.per_target_cov \
    --bed ${bed} \
    --out ${outdir}/not_merged/${prefix}.QCstats.tsv


# pileup-bcftools

#samtools mpileup -t DP -t SP -t AD --VCF --uncompressed --max-depth 200000 \
#    -o ${outdir}/not_merged/${prefix}.mpileup.vcf \
#    -f ${refFa} -l ${bed} --min-MQ 30 --min-BQ 20 \
#    ${outdir}/not_merged/${prefix}.sorted.bam
#bcftools call -v -m -P 0.00001 -T ${bed} -A ${outdir}/not_merged/${prefix}.mpileup.vcf > ${outdir}/not_merged/${prefix}.mpileup.filter.vcf

# vardict

${VardictBin}/VarDict -G ${refFa} \
    -N "SAMPLE" \
    -b "${outdir}/not_merged/${prefix}.sorted.bam" \
    -th 4 -O 3 -Q 3 -F 0 -p -a 20:0.2 \
    --nosv -f 0.00001 -c 1 -S 2 -E 3 -g 4 ${outdir}/mutations/${bed_name}.amplicon.bed | \
    ${VardictBin}/teststrandbias.R | \
    ${VardictBin}/var2vcf_valid.pl \
        -A -N "SAMPLE" -d 1000 -v 3 -f 0.00001 -P 0  > ${outdir}/not_merged/${prefix}.vardict.vcf

# filter

bgzip ${outdir}/not_merged/${prefix}.vardict.vcf && tabix ${outdir}/not_merged/${prefix}.vardict.vcf.gz
bcftools view -i "FORMAT/VD>=10 && %FILTER='PASS'" ${outdir}/not_merged/${prefix}.vardict.vcf.gz > ${outdir}/not_merged/${prefix}.vardict.filter.vcf

#bedtools intersect -wb -a ${outdir}/not_merged/${prefix}.vardict.filter.vcf -b ${site} | awk '$4==$14&&$5==$15' > ${outdir}/not_merged/${prefix}.vardict.filter.intersect.tsv

#awk '{if($1~/^#/){print $0} else{split($10,a,":"); if(($4=="G"&&$5=="A"&&a[5]>=0.0003)||($4=="C"&&$5=="T"&&a[5]>=0.0003)){print $0}}}' ${outdir}/not_merged/${prefix}.vardict.filter.intersect.tsv > ${outdir}/not_merged/${prefix}.vardict.filter.intersect.filter_transition.tsv


