#!/usr/bin/bash
# -*- coding: utf-8 -*-
# @Author  : zhangsiwen


# <<<--- variables --->>>

shellname=$0
scriptDir=`dirname ${shellname}`
. ${scriptDir}/src/env.sh
#refFa='/mnt/user/errand/pipeline/MC3/ref/Homo_sapiens_assembly19.fasta'
refFa='/mnt/ddngs/zhangsw/database/GRCh37_clinical/nochrY/Homo_sapiens_assembly19.nochrY.fasta'
VardictBin=/mnt/ddngs/zhangsw/software/VarDict-1.8.3/bin
fastp=/mnt/ddngs/zhangsw/bin/fastp
usearch=/mnt/user/share/bin/usearch
picard=/mnt/ddngs/zhangsw/miniconda3/bin/picard
Rscript='/mnt/ddngs/zhangsw/miniconda3/bin/Rscript'
bcftools=/mnt/ddngs/zhangsw/miniconda3/envs/bcftools/bin/bcftools
bedtools=/mnt/user/share/bin/bedtools
outdir=`pwd`'/outdir'
prefix='sample'
ad=1
af=0.00001
thread=4
clean=low
force=false


# usage
function usage {
    printf -- "Usage:\n    sh ${shellname} -f R1.fq -r R2.fq -l bed -s site.bed -v refFa -o outdir -p prefix -n 3 -t thread -c low \n\n";
    printf -- "Required arguments: \n";
    printf -- "    -f|--forward       R1.fq \n";
    printf -- "    -r|--reverse       R2.fq \n";
    printf -- "    -l|--bed           bed region \n";
    printf -- "    -s|--site          site coordinate bed \n";
    printf -- "Optional arguments: \n";
    printf -- "    -v|--refFa         reference fasta, default '${refFa}' \n"
    printf -- "    -o|--outdir        output path, default '${outdir}' \n";
    printf -- "    -p|--prefix        output prefix, default '${prefix}' \n";
    printf -- "    -n|--ad            minimum allele depth to be reported, default '${ad}' \n";
    printf -- "    -a|--af            minimum allele frequency to be reported, default '${af}' \n";
    printf -- "    -t|--thread        thread, default '${thread}' \n";
    printf -- "    -c|--clean         if clean intermediate files, levels: none|low|high, default '${clean}' \n";
    printf -- "                       none: keep all results,\n";
    printf -- "                       low: remove fastp/*fq, fastq_merge/*merged.fq, mapping/*sam, mapping/unsorted.bam, \n";
    printf -- "                       high: remove fastp/, fastq_merge/, mapping/, \n";
    printf -- "    -f|--force         force to over-write output \n";
    printf -- "    -h|--help          help \n";
    printf -- "@Author: Siwen Zhang \n\n";
    exit 0
}


# Arguments passing
TEMP=`getopt -o f:r:l:s:v:o:p:n:a:t:c:fh --long forward:,reverse:,bed:,site:,refFa:,outdir:,prefix:,ad:,af:,thread:,clean:,force,help \
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
        -n|--ad)
            ad="$2"; shift 2 ;;
        -a|--af)
            af="$2"; shift 2 ;;
        -t|--thread)
            thread="$2"; shift 2 ;;
        -c|--clean)
            clean="$2"; shift 2;;
        -f|--force)
            force=true ; shift ;;
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


# overwrite results
checkResults () {
    if $1 
    then true
    else if [ -f $2 ]
        then false
        else true
        fi
    fi
}


# arguments not null
echo "R1: ${forward}"
echo "R2: ${reverse}"
echo "BED: ${bed}"
echo "SITE: ${site}"
if [ -f ${forward} ] && [ -f ${reverse} ] && [ -f ${bed} ] && [ -f ${site} ]
then echo "Files found."
else echo "Files NOT found."
    usage
fi


# pre-process

mkdir -p ${outdir}/fastp

if checkResults ${force} ${outdir}/fastp/${prefix}.read_stat.txt
then
${fastp} --in1 ${forward} \
    --in2 ${reverse} \
    --out1 ${outdir}/fastp/${prefix}_R1.fq \
    --out2 ${outdir}/fastp/${prefix}_R2.fq \
    --json ${outdir}/fastp/${prefix}.json \
    --html ${outdir}/fastp/${prefix}.html \
    --adapter_fasta ${scriptDir}/src/adapters.fasta \
    --detect_adapter_for_pe \
    --average_qual 20 --cut_mean_quality 20 \
    --thread ${thread} --reads_to_process 12000000 > ${outdir}/fastp/log.fastp 2>&1

python ${scriptDir}/src/read_stat.v1.py \
    --fastp ${outdir}/fastp/${prefix}.json \
    --out ${outdir}/fastp/${prefix}.read_stat.txt
fi


# pear merge fq1 and fq2

mkdir -p ${outdir}/fastq_merge 

# usearch merge fq1 and fq2
#https://drive5.com/usearch/manual/cmds_all.html
#https://drive5.com/usearch/manual/exp_errs.html

# debug mode add:  -alnout ${outdir}/fastq_merge/${prefix}.merged.aln 
# use clean mode to decide if use debug mode, 20230317

if [ ${clean} == "none" ]
then usearch_aln="-alnout ${outdir}/fastq_merge/${prefix}.merged.aln"
    usearch_notmerged="-fastqout_notmerged_fwd ${outdir}/fastq_merge/${prefix}.notmerged_R1.fq -fastqout_notmerged_rev ${outdir}/fastq_merge/${prefix}.notmerged_R2.fq "
else usearch_aln=""
    usearch_notmerged=""
fi

if checkResults ${force} ${outdir}/fastq_merge/${prefix}.merged.fq
then
${usearch} -fastq_mergepairs ${outdir}/fastp/${prefix}_R1.fq \
    -reverse ${outdir}/fastp/${prefix}_R2.fq \
    -fastqout ${outdir}/fastq_merge/${prefix}.merged.fq \
    -tabbedout ${outdir}/fastq_merge/${prefix}.merged.tsv \
    -threads ${thread} \
    -fastq_maxdiffs 10 -fastq_pctid 90 -fastq_merge_maxee 0.01 \
    -fastq_minovlen 20 -fastq_minmergelen 20 -fastq_maxmergelen 200 -fastq_minlen 20 \
    ${usearch_aln}  ${usearch_notmerged}  2> ${outdir}/fastq_merge/${prefix}.merged.stats
fi

sed -n '10,11p' ${outdir}/fastq_merge/${prefix}.merged.stats|awk 'NR==1{total=$1;print "Pair_Reads\t"$1}NR==2{n=$1;print "Merged_Reads\t"$1}END{print "Merged_Reads_Pct\t"n/total}' > ${outdir}/fastq_merge/${prefix}.merged.stats.tsv


if [ -f ${outdir}/fastq_merge/R1.fq ]
then rm ${outdir}/fastq_merge/R1.fq ${outdir}/fastq_merge/R2.fq
fi


# ----- merged reads ----- #

# mapping

mkdir -p ${outdir}/mapping

if checkResults ${force} ${outdir}/mapping/${prefix}.sorted.bam
then
bwa mem -t ${thread} -M \
    -R "@RG\\tID:NULL\\tSM:${prefix}\\tLB:NULL\\tPU:NULL\\tPL:illumina\\tCN:prec_sci" \
    ${refFa} ${outdir}/fastq_merge/${prefix}.merged.fq 2> log.txt 1> ${outdir}/mapping/${prefix}.sam

samtools view -Sb -@ ${thread} -t ${refFa}.fai ${outdir}/mapping/${prefix}.sam > ${outdir}/mapping/${prefix}.bam

sambamba sort -t ${thread} -l 1 -o ${outdir}/mapping/${prefix}.sorted.bam ${outdir}/mapping/${prefix}.bam

samtools index -@ ${thread} ${outdir}/mapping/${prefix}.sorted.bam ${outdir}/mapping/${prefix}.sorted.bam.bai 
fi

# stats

mkdir -p ${outdir}/stats

# count base number for each region
#https://broadinstitute.github.io/picard/

if [ -f ${refFa}.dict ]
then refDict=${refFa}.dict
else refDict=`echo ${refFa}|sed 's/fasta$//'|sed 's/fa$//'`"dict"
fi

bed_name=`basename -s '.bed' ${bed}`

if checkResults ${force} ${outdir}/stats/${prefix}.QCstats.tsv
then
${picard} BedToIntervalList \
    -I ${bed} \
    -O ${outdir}/stats/${bed_name}.intervalList \
    -SD ${refDict}

${picard} CollectTargetedPcrMetrics \
    I=${outdir}/mapping/${prefix}.sorted.bam \
    O=${outdir}/stats/${prefix}.CollectTargetedPcrMetrics \
    R=${refFa} \
    AMPLICON_INTERVALS=${outdir}/stats/${bed_name}.intervalList \
    TARGET_INTERVALS=${outdir}/stats/${bed_name}.intervalList \
    PER_TARGET_COVERAGE=${outdir}/stats/${prefix}.per_target_cov \
    PER_BASE_COVERAGE=${outdir}/stats/${prefix}.per_base_cov \
    MINIMUM_BASE_QUALITY=20 \
    NEAR_DISTANCE=100 \
    COVERAGE_CAP=100
    MINIMUM_MAPPING_QUALITY=30 

python ${scriptDir}/src/mapping_qc_stats.py \
    --metrics ${outdir}/stats/${prefix}.CollectTargetedPcrMetrics \
    --per_base_cov ${outdir}/stats/${prefix}.per_base_cov \
    --per_target_cov ${outdir}/stats/${prefix}.per_target_cov \
    --bed ${bed} \
    --out ${outdir}/stats/${prefix}.QCstats.tsv
fi

${Rscript} --vanilla ${scriptDir}/src/boxplot_with_quantiles.r \
    --in ${outdir}/stats/${prefix}.per_target_cov \
    --out ${prefix}.per_target_cov \
    --dir ${outdir}/stats/ \
    --key mean_coverage

cat ${outdir}/fastp/${prefix}.read_stat.txt ${outdir}/fastq_merge/${prefix}.merged.stats.tsv ${outdir}/stats/${prefix}.QCstats.tsv > ${outdir}/${prefix}.all.QCstats.tsv



# mutations

mkdir ${outdir}/mutations

# vardict

awk 'BEGIN{FS="\t";OFS="\t"}{if($1!=""){print $1,$2-20,$3+20,$4,"0",".",$2,$3}}' ${bed} > ${outdir}/mutations/${bed_name}.amplicon.bed

if checkResults ${force} ${outdir}/mutations/${prefix}.vardict.vcf
then
${VardictBin}/VarDict -G ${refFa} \
    -N "SAMPLE" \
    -b "${outdir}/mapping/${prefix}.sorted.bam" \
    -th ${thread} -O 20 -Q 20 -F 0 -p -a 20:0.2 -q 25 -o 3 \
    --nosv -f ${af} -r 1 -c 1 -S 2 -E 3 -g 4 ${outdir}/mutations/${bed_name}.amplicon.bed | \
    ${VardictBin}/teststrandbias.R | \
    ${VardictBin}/var2vcf_valid.pl \
        -A -N "SAMPLE" -d 1000 -v ${ad} -f ${af} -P 0 > ${outdir}/mutations/${prefix}.vardict.vcf
fi


# filter

if checkResults ${force} ${outdir}/mutations/${prefix}.vardict.filter.vcf
then
bgzip -f ${outdir}/mutations/${prefix}.vardict.vcf && tabix ${outdir}/mutations/${prefix}.vardict.vcf.gz
${bcftools} norm -a -O v -m - ${outdir}/mutations/${prefix}.vardict.vcf.gz| ${bedtools} intersect -header -a - -b ${bed} |uniq > ${outdir}/mutations/${prefix}.vardict.filter.vcf
fi

perl ${scriptDir}/src/vcf2maf.pl \
    --input-vcf ${outdir}/mutations/${prefix}.vardict.filter.vcf \
    --output-maf ${outdir}/mutations/${prefix}.vardict.filter.maf \
    --tumor-id SAMPLE \
    --inhibit-vep --ref-fasta ${refFa} \
    --filter-vcf 0 --ncbi-build GRCh37 --retain-info REF,ALT,COSMIC,CENTERS,CONTEXT,DBVS


if checkResults ${force} ${outdir}/mutations/${prefix}.vardict.filter.intersect.result.tsv
then
awk -v ad=${ad} 'BEGIN{FS="\t";OFS="\t"}
    NR==FNR{if($1!~"#"&&$1!="Hugo_Symbol"){if($42>=ad){af=$42/$40}else{af=0}; if(a[$5":"$6":"$11":"$13]){if(freq[$5":"$6":"$11":"$13]<af){a[$5":"$6":"$11":"$13]=$11">"$13"\t"$40"\t"$42"\t"af; freq[$5":"$6":"$11":"$13]=af}} else{a[$5":"$6":"$11":"$13]=$11">"$13"\t"$40"\t"$42"\t"af; freq[$5":"$6":"$11":"$13]=af}}}
    NR>FNR{if($1=="#Chromosome"){print $0,"detected","depth","reads","vaf"} else{if(a[$1":"$2+1":"$4":"$5]){print $0,a[$1":"$2+1":"$4":"$5]}else{print $0,$4">"$5,"0","0","0"}}}' ${outdir}/mutations/${prefix}.vardict.filter.maf ${site}| \
    awk 'BEGIN{FS="\t";OFS="\t"}
    NR==FNR{a[$1":"$2]=$4}
    NR>FNR{if($10==0){print $1,$2,$3,$4,$5,$6,$7,$8,$4">"$5,a[$1":"$2+1],"0","0"} else{print $0}}' ${outdir}/stats/${prefix}.per_base_cov - > ${outdir}/mutations/${prefix}.vardict.filter.intersect.tsv

awk 'BEGIN{FS="\t";OFS="\t"}NR==1{print $0,"result"}NR>1{if($NF>=0.0003){print $0,"Positive"}else{print $0,"Negative"}}' ${outdir}/mutations/${prefix}.vardict.filter.intersect.tsv > ${outdir}/mutations/${prefix}.vardict.filter.intersect.result.tsv
fi


# final clean intermediate files

echo "CLEAN mode: "${clean}
if [ ${clean} == "high" ] 
then echo "Cleaning: ${outdir}/fastp/  fastq_merge/  ${outdir}/mapping/"
    rm -r ${outdir}/fastp/  fastq_merge/  ${outdir}/mapping/
elif [ ${clean} == "low" ] 
then echo "Cleaning: ${outdir}/fastp/*fq  ${outdir}/fastq_merge/*merged.fq  ${outdir}/mapping/${prefix}.bam  ${outdir}/mapping/${prefix}.sam"
    rm ${outdir}/fastp/*fq  ${outdir}/fastq_merge/*merged.fq  ${outdir}/mapping/${prefix}.bam  ${outdir}/mapping/${prefix}.sam
fi

