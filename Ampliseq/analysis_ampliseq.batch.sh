#!/usr/bin/bash
# -*- coding: utf-8 -*-
# @Author  : zhangsiwen


# <<<--- variables --->>>

shellname=$0
scriptDir=`dirname ${shellname}`
. ${scriptDir}/src/env.sh
#refFa='/mnt/user/errand/pipeline/MC3/ref/Homo_sapiens_assembly19.fasta'
refFa='/mnt/ddngs/zhangsw/database/GRCh37_clinical/nochrY/Homo_sapiens_assembly19.nochrY.fasta'
Rscript='/mnt/ddngs/zhangsw/miniconda3/bin/Rscript'
outdir=`pwd`'/Outdir'
thread=4
clean=low
qhost=clinical2


# usage
function usage {
    printf -- "Usage:\n    bash ${shellname} -i config -l bed -s site.bed -v refFa -o Outdir -p prefix -t thread -c low -q clinical.q \n\n";
    printf -- "Required arguments: \n";
    printf -- "    -i|--input         input config file, tab-delimited, Group column only accept 'control' or 'test',
                       Sample column will be used as output dir name and file prefix name\n";
    printf -- "    -l|--bed           bed region \n";
    printf -- "    -s|--site          site coordinate bed \n";
    printf -- "Optional arguments: \n";
    printf -- "    -v|--refFa         reference fasta, default '${refFa}' \n"
    printf -- "    -o|--outdir        output path, default '${outdir}' \n";
    printf -- "    -t|--thread        thread, default '${thread}' \n";
    printf -- "    -c|--clean         if clean intermediate files, levels: none|low|high, default '${clean}' \n";
    printf -- "                       none: keep all results,\n";
    printf -- "                       low: remove fastp/*fq, fastq_merge/*merged.fq, mapping/*sam, mapping/unsorted.bam, \n";
    printf -- "                       high: remove fastp/, fastq_merge/, mapping/, \n";
    printf -- "    -q|--qhost         sge host name, default '${qhost}'  \n\n";
    printf -- "    -h|--help          help \n\n";
    printf -- "Input file example (Tab-delimited):
    Sample	R1	R2	Group
    N1	N1_R1.fq	N1_R2.fq	control
    N2	N2_R1.fq	N2_R2.fq	control
    N3	N3_R2.fq	N3_R3.fq	control
    N4	N4_R2.fq	N4_R3.fq	control
    N5	N5_R3.fq	N5_R4.fq	control
    S1	S1_R1.fq	S1_R2.fq	test
    S2	S2_R1.fq	S2_R2.fq	test\n\n"
    printf -- "@Author: Siwen Zhang \n\n";
    exit 0
}


# Arguments passing
TEMP=`getopt -o i:l:s:v:o:p:t:c:q:h --long input:,bed:,site:,refFa:,outdir:,prefix:,thread:,clean:,qhost:,help \
     -n ${shellname} -- "$@"`

# no arguments detected
if [ $? != 0 ] || [ $# == 0 ]
    then usage
fi

# set arguments
eval set -- "$TEMP"
while true ; do
    case "$1" in
        -i|--input)
            input="$2"; shift 2 ;;
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
        -c|--clean)
            clean="$2"; shift 2;;
        -q|--qhost)
            qhost="$2"; shift 2;;
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


# a function to check input files

currentDir=`pwd`

function checkFile () {
    if [[ $1 == /* ]]
    then if [ -f $1 ]
        then echo $1
        else echo "File not found: "$1; exit 0;
        fi
    else if [ -f ${currentDir}"/"$1 ]
        then echo ${currentDir}"/"$1
        else echo "File not found: "$1; exit 0;
        fi
    fi
}

function checkDir () {
    if [[ $1 == /* ]]
    then echo $1
    else echo ${currentDir}"/"$1
    fi
}


# check sample input

Input=`checkFile ${input}`
Bed=`checkFile ${bed}`
Site=`checkFile ${site}`
Outdir=`checkDir ${outdir}`
mkdir -p ${Outdir}/run_bash
echo "Input: ${Input}"
echo "BED: ${Bed}"
echo "SITE: ${Site}"
echo "Outdir: ${Outdir}"
echo "CLEAN mode: "${clean}
d=`date +"%Y_%m_%d_%H_%M_%S"`

if [ -f "${Input}" ] && [ -f "${Bed}" ] && [ -f "${Site}" ]
then echo "Files found."
    echo "Total sample number: "`awk '$1!="Sample"&&$1!=""&&$4!="Group"' ${Input}|wc -l`
    awk '$4=="control"||$4=="Control"' ${Input} > ${Outdir}/sample_list.control.${d}
    control_sample=`cat ${Outdir}/sample_list.control.${d}|wc -l`
    echo "Control sample number: "${control_sample}
    awk '$4=="test"||$4=="Test"' ${Input} > ${Outdir}/sample_list.test.${d}
    test_sample=`cat ${Outdir}/sample_list.test.${d}|wc -l`
    echo "Test sample number: "${test_sample}
else echo "Files NOT found."
    usage
fi


# run bash directory
cd ${Outdir}/run_bash

# job for each sample

awk '$1!="Sample"&&$1!=""&&$4!="Group"' ${Input}|\
while read sample r1 r2 group
do 
R1=`checkFile ${r1}`
R2=`checkFile ${r2}`
echo "${sample} R1: ${R1}"
echo "${sample} R2: ${R2}"
if [ -f "${R1}" ] && [ -f "${R2}" ]
then echo "${sample} Fastq files found."
else echo "${sample} Fastq files NOT found."
    exit 0;
fi
if [ "${group}" == "Test" ] || [ "${group}" == "test" ] 
then ad=10; af=0.000001
else ad=1; af=0.000001
fi
echo "${sample} Set AD: "${ad}" and AF: "${af}
echo ". ${scriptDir}/src/env.sh" > ${Outdir}/run_bash/run.${sample}.sh
echo "bash ${scriptDir}/analysis_pipeline.sh -f ${R1} -r ${R2} -l ${Bed} -s ${Site} -v ${refFa} -o ${Outdir}/${group}-${sample} -p ${sample} -n ${ad} -a ${af} -t ${thread} -c ${clean} " >> ${Outdir}/run_bash/run.${sample}.sh
qsub -cwd -l mf=12g,num_proc=${thread} -q ${qhost} -N "J"${sample} ${Outdir}/run_bash/run.${sample}.sh
done


# polishing for all control and test samples
mkdir -p ${Outdir}/polishing_${d}/

# qsub jobid
jobid=$(awk 'BEGIN{ORS=","}$1!="Sample"&&$1!=""&&$4!="Group"{print "J"$1}' ${Input}|sed -e 's/,$//') 
echo "QSUB jobid: ${jobid}"

# polishing run bash
echo ". ${scriptDir}/src/env.sh" > ${Outdir}/run_bash/run.polishing_${d}.sh
echo "
for i in \`ls ${Outdir}/*/mutations/*filter.intersect.result.tsv\`; 
do 
awk -v sample=\`basename -s .vardict.filter.intersect.result.tsv \${i}\` 'NR==1{print \"mutation_id\\t\"sample}NR>1{print \$1\":\"\$2+1\":\"\$4\":\"\$5\"\\t\"\$10}' \${i} > \${i}.DP; 
awk -v sample=\`basename -s .vardict.filter.intersect.result.tsv \${i}\` 'NR==1{print \"mutation_id\\t\"sample}NR>1{print \$1\":\"\$2+1\":\"\$4\":\"\$5\"\\t\"\$12}' \${i} > \${i}.VAF; 
done
" >> ${Outdir}/run_bash/run.polishing_${d}.sh
echo "control_sample=\`ls -d ${Outdir}/[Cc]ontrol-*/|wc -l\`
test_sample=\`ls -d ${Outdir}/[Tt]est-*/|wc -l\`
all_sample=\`ls -d ${Outdir}/*/|wc -l\`
paste ${Outdir}/[Cc]ontrol-*/mutations/*.vardict.filter.intersect.result.tsv.DP|cut -f\`bash ${scriptDir}/src/get_cut_field.sh 2 \${control_sample} start\` > ${Outdir}/polishing_${d}/allsample.control.DP.tsv
paste ${Outdir}/[Cc]ontrol-*/mutations/*.vardict.filter.intersect.result.tsv.VAF|cut -f\`bash ${scriptDir}/src/get_cut_field.sh 2 \${control_sample} start\` > ${Outdir}/polishing_${d}/allsample.control.VAF.tsv
paste ${Outdir}/[Tt]est-*/mutations/*.vardict.filter.intersect.result.tsv.DP|cut -f\`bash ${scriptDir}/src/get_cut_field.sh 2 \${test_sample} start\` > ${Outdir}/polishing_${d}/allsample.test.DP.tsv
paste ${Outdir}/[Tt]est-*/mutations/*.vardict.filter.intersect.result.tsv.VAF|cut -f\`bash ${scriptDir}/src/get_cut_field.sh 2 \${test_sample} start\` > ${Outdir}/polishing_${d}/allsample.test.VAF.tsv
paste ${Outdir}/*/mutations/*.vardict.filter.intersect.result.tsv.DP|cut -f\`bash ${scriptDir}/src/get_cut_field.sh 2 \${all_sample} start\` > ${Outdir}/polishing_${d}/allsample.DP.tsv
${Rscript} --vanilla ${scriptDir}/src/boxplot_with_quantiles.r --input ${Outdir}/polishing_${d}/allsample.DP.tsv --out allsample.DP --dir ${Outdir}/polishing_${d}/ 
" >> ${Outdir}/run_bash/run.polishing_${d}.sh

echo "${Rscript} --vanilla ${scriptDir}/src/polishing_rebuild.r --nctvaf ${Outdir}/polishing_${d}/allsample.control.VAF.tsv --nctdp ${Outdir}/polishing_${d}/allsample.control.DP.tsv --testvaf ${Outdir}/polishing_${d}/allsample.test.VAF.tsv --testdp ${Outdir}/polishing_${d}/allsample.test.DP.tsv --plot all --out ${Outdir}/polishing_${d}/all_sample.final --pval 0.05
paste ${Outdir}/polishing_${d}/all_sample.final.polishing.tsv ${Outdir}/polishing_${d}/all_sample.final.model_fit.tsv > ${Outdir}/polishing_${d}/all_sample.final.result.tsv" >> ${Outdir}/run_bash/run.polishing_${d}.sh

qsub -cwd -l mf=12g -q ${qhost} -hold_jid ${jobid} -N MRDpolishing ${Outdir}/run_bash/run.polishing_${d}.sh
