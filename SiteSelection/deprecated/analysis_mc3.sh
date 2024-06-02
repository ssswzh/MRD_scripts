#!/usr/bin/bash
# path
# /public2/home/lizhe/zhangsw/jinyun/project/PRJZ_001_Pan_gynecologic_wes_mutation_research/Standard_Analysis/Chemotherapy/
# /public2/home/lizhe/zhangsw/jinyun/project/PRJZ_001_Pan_gynecologic_wes_mutation_research/Standard_Analysis/PAPRI/20210830_WES

# jupiter
# /mnt/ddngs/zhangshouwei/project/PRJZ_001_Pan_gynecologic_wes_mutation_research/PAPRI/Basis_results/Results/MC3

# 雄衡2
# /public/home/zhangsw01/project/zhangsw/jinyun/project/PRJZ_001_Pan_gynecologic_wes_mutation_research/Standard_Analysis/Chemotherapy


# <<<--- variables --->>>

#tumorBam='bams/P04_T.raw.sorted.markdup.realn.recal.bam'
#tumorSampleName='P04_T'
#somaticVcf='snv_cnv/P04.merge_vcfs.final.flt.vep.vcf'
#germlineVcf='snv_cnv/P04_N.germline.final.vep.vcf'
#cnv='snv_cnv/P04_T.raw.sorted.markdup.realn.call.cns'
#out='P04'
#outdir='/mnt/ddngs/zhangsw/project/MRD/OV_WES/P04/'
scriptDir='/mnt/ddngs/zhangsw/project/MRD/OV_WES/scripts'
refFa='/mnt/user/errand/pipeline/MC3/ref/Homo_sapiens_assembly19.fasta'
bed='/mnt/ddngs/zhangsw/project/MRD/OV_WES/bed/Agilent_v6r2.Covered.bed'
interval='/mnt/user/errand/pipeline/MC3/ref/Agilent.SureSelect_All_Exon-v6_r2.interval_list'
annoDir='/mnt/ddngs/zhangsw/database/GRCh37'
snp_pileup='/mnt/user/errand/pipeline/OncoDeficiencyPro/tools/snppileup/htstools-master/snp-pileup'
#run_facet='/mnt/ddngs/lvr/script/HRD/OncoDeficiencyPro/script/run_facet.lr.R'
#facets_score='/mnt/ddngs/lvr/script/HRD/OncoDeficiencyPro/script/get_scarHRD_score_plus.R'
Rscript='/mnt/ddngs/lvr/Miniconda3/miniconda3/envs/py2.7/bin/Rscript'


# Arguments passing
shellname=$0
TEMP=`getopt -o b:n:s:v:g:c:o:p:h --long tbam:,nbam:,tsample:,somaticVcf:,germlineVcf:,cnv:,outdir:,prefix:,help \
     -n ${shellname} -- "$@"`

# no arguments detected
if [ $? != 0 ] || [ $# == 0 ]
    then
        printf "Error: please try again with correct arguments.\n\n";
        printf "Usage:\nsh $0  \n\n" >&2 ;
        printf "To see detailed arguments, use '-h' or '--help' \n\n";
        exit 1 ;
fi

# set arguments
eval set -- "$TEMP"
while true ; do
    case "$1" in
        -b|--tbam)
            tbam="$2"; shift 2;;
        -n|--nbam)
            nbam="$2"; shift 2;;
        -s|--tsample)
            tsample="$2"; shift 2;;
        -v|--somaticVcf)
            somaticVcf="$2"; shift 2;;
        -g|--germlineVcf)
            germlineVcf="$2"; shift 2 ;;
        -c|--cnv)
            cnv="$2"; shift 2 ;;
        -o|--outdir)
            outdir="$2"; shift 2 ;;
        -p|--prefix)
            prefix="$2"; shift 2 ;;
        -h|--help)
            printf -- "Usage:\nsh ${shellname} -b tumorBam -n normalBam -s tumorSampleName -v somaticVcf -g germlineVcf -c cnv -o outdir -p prefix\n\n";
            printf -- "-b|--tbam          tumor bam file \n";
            printf -- "-n|--nbam          normal bam file \n";
            printf -- "-s|--tsample       tumor sample name \n";
            printf -- "-v|--somaticVcf    tumor somatic vcf \n"
            printf -- "-g|--germlineVcf   normal germline vcf \n";
            printf -- "-c|--cnv           cnv result \n";
            printf -- "-o|--outdir        output path \n";
            printf -- "-p|--prefix        output prefix \n";
            printf -- "-h|--help          help \n";
            printf -- "@Author: Siwen Zhang, zhang.siwen@puruijizhun.com \n\n";
            exit 0;; #shift ; break; exit 1;;
        --)
            shift ; break; exit 1;;
        ?)
            printf -- 'Usage: \n';
            printf -- 'sh ${shellname} -b tumorBam -n normalBam -s tumorSampleName -v somaticVcf -g germlineVcf -c cnv -o outdir -p prefix \n\n';
            printf -- '@Author: Siwen Zhang, zhang.siwen@puruijizhun.com \n\n';
            exit 0;; #shift ; break; exit 1;;
        *)
            printf -- 'Usage: ';
            printf -- 'sh ${shellname} -b tumorBam -n normalBam -s tumorSampleName -v somaticVcf -g germlineVcf -c cnv -o outdir -p prefix \n\n'
            printf -- '@Author: Siwen Zhang, zhang.siwen@puruijizhun.com \n\n';
            exit 0;; #exit 1;;
    esac
done



# <<<--- SNV indel calling --->>>

## filter blacklist (Duke, DAC, repeatMasker, simpleRepeat)

### somatic

mkdir -p ${outdir}/snv_cnv/
bedtools subtract -nonamecheck -header -a ${somaticVcf} -b ${annoDir}/wgEncodeDacMapabilityConsensusExcludable.bed.gz |\
    bedtools subtract -nonamecheck -header -a - -b ${annoDir}/wgEncodeDukeMapabilityRegionsExcludable.bed.gz |\
    bedtools subtract -nonamecheck -header -a - -b ${annoDir}/hg19_repeatMasker.bed > ${outdir}/snv_cnv/${prefix}.somatic.filterBlacklist.vcf

bcftools norm -O v -o ${outdir}/snv_cnv/${prefix}.somatic.filterBlacklist.norm.vcf -m + ${outdir}/snv_cnv/${prefix}.somatic.filterBlacklist.vcf

bcftools view --types snps ${outdir}/snv_cnv/${prefix}.somatic.filterBlacklist.norm.vcf > ${outdir}/snv_cnv/${prefix}.somatic.filterBlacklist.norm.onlysnp.vcf 

# germline

bedtools subtract -nonamecheck -header -a ${germlineVcf} -b ${annoDir}/wgEncodeDacMapabilityConsensusExcludable.bed.gz |\
    bedtools subtract -nonamecheck -header -a - -b ${annoDir}/wgEncodeDukeMapabilityRegionsExcludable.bed.gz |\
    bedtools subtract -nonamecheck -header -a - -b ${annoDir}/hg19_repeatMasker.bed > ${outdir}/snv_cnv/${prefix}.germline.filterBlacklist.vcf

bcftools norm -O v -o ${outdir}/snv_cnv/${prefix}.germline.filterBlacklist.norm.vcf -m + ${outdir}/snv_cnv/${prefix}.germline.filterBlacklist.vcf




# <<<--- Copy Number Variation --->>>

## germline variants, calculate BAF

#python ${scriptDir}/calculate_baf.py \
#    --vcf ${outdir}/snv_cnv/${prefix}.germline.filterBlacklist.norm.vcf \
#    --out ${outdir}/snv_cnv/${prefix}.germline.filterBlacklist.norm.baf \
#    --fmt tsv --source other

## extract somatic BAF from germline tsv

#/mnt/ddngs/zhangsw/miniconda3/envs/bam-readcount/bin/bam-readcount \
#    --min-mapping-quality 1 --min-base-quality 20 --max-warnings 0 \
#    -f ${refFa} -l ${bed} \
#    ${tbam} > ${outdir}/${prefix}.somatic.bam-readcount

#python ${scriptDir}/readcount_extraction.py \
#    --vcf ${outdir}/snv_cnv/${prefix}.somatic.filterBlacklist.norm.onlysnp.vcf \
#    --count ${outdir}/${prefix}.somatic.bam-readcount \
#    --out ${outdir}/snv_cnv/${prefix}.tumor.baf.tsv \
#    --mode vaf --soft bam-readcount 


## run facets

mkdir ${outdir}/facets
${snp_pileup} -g -q1 -Q20 ${outdir}/snv_cnv/${prefix}.germline.filterBlacklist.norm.vcf ${outdir}/facets/${prefix}.snp-pileup ${nbam} ${tbam}

${Rscript} ${scriptDir}/run_facet.lr.zsw.R ${outdir}/facets/${prefix}.snp-pileup.gz ${prefix} ${outdir}/facets/${prefix}






## extract logR from BAF file

#python ${scriptDir}/extract_mutation_cnv.py \
#    --baf ${outdir}/ascat/${prefix}.germline.filterBlacklist.norm.baf.tsv \
#    --cnv ${cnv} \
#    --out ${outdir}/ascat/${prefix}.tumor.logR.tsv \
#    --source cnvkit


## normal logR, calculated by dividing median depth

#awk 'BEGIN{FS="\t";OFS="\t"}NR==1{print $0}NR>1{print $1,$2,$3,0}' ${outdir}/ascat/${prefix}.tumor.logR.tsv > ${outdir}/ascat/${prefix}.normal.logR.tsv

#picard CollectHsMetrics \
#    I=${outdir}/bams/${prefix}_N.raw.sorted.markdup.realn.recal.bam \
#    O=${outdir}/depth/${prefix}.normal.coverage \
#    PER_TARGET_COVERAGE=${outdir}/depth/${prefix}.normal.perTarget.coverage \
#    PER_BASE_COVERAGE=${outdir}/depth/${prefix}.normal.perBase.depth \
#    BAIT_INTERVALS=${interval} \
#    TARGET_INTERVALS=${interval} \
#    R=${refFa} 

#python ${scriptDir}/extract_mutation_cnv.py \
#    --baf ${outdir}/ascat/${prefix}.germline.filterBlacklist.norm.baf.tsv \
#    --cnv ${outdir}/depth/${prefix}.normal.perTarget.coverage \
#    --out ${outdir}/ascat/${prefix}.tumor.logR.tsv \
#    --source picard



## run ASCAT

#Rscript --vanilla ${scriptDir}/ascat_cnv.r \
#    --somaticcnv ${outdir}/ascat/${prefix}.tumor.logR.tsv \
#    --somaticbaf ${outdir}/ascat/${prefix}.tumor.baf.tsv \
#    --germlinecnv ${outdir}/ascat/${prefix}.normal.logR.tsv \
#    --germlinebaf ${outdir}/ascat/${prefix}.germline.filterBlacklist.norm.baf.tsv \
#    --outdir ${outdir}/ascat/ \
#    --out ${prefix}.cnv.ascat



# <<<--- PyClone cluster and phylogenetic tree construction --->>>

## prepare data

mkdir -p ${outdir}/pyclone_analysis/

perl ${scriptDir}/vcf2maf.pl \
    --input-vcf ${outdir}/snv_cnv/${prefix}.somatic.filterBlacklist.norm.vcf \
    --output-maf ${outdir}/snv_cnv/${prefix}.somatic.filterBlacklist.norm.maf \
    --tumor-id PRIMARY --normal-id NORMAL \
    --custom-enst ${scriptDir}/GCF_000001405.25_GRCh37.p13_genomic.selected_transcript.txt \
    --inhibit-vep --ref-fasta ${refFa} \
    --filter-vcf 0 --ncbi-build GRCh37 --retain-info COSMIC,CENTERS,CONTEXT,DBVS

python ${scriptDir}/pyclone_dat_generation.py \
    --infile ${outdir}/snv_cnv/${prefix}.somatic.filterBlacklist.norm.maf \
    --cnv ${outdir}/facets/${prefix}/${prefix}.4.segment.final.txt \
    --out ${outdir}/pyclone_analysis/${prefix}.pyclone_dat.tsv \
    --gender female --sampleid ${prefix} --depth 0 --ad 0 --af 0 --source facets

## run new pyclone

/mnt/ddngs/zhangsw/miniconda3/envs/pyclone-vi/bin/pyclone-vi fit \
    -i ${outdir}/pyclone_analysis/${prefix}.pyclone_dat.tsv \
    -o ${outdir}/pyclone_analysis/${prefix}.pyclone_dat.h5 \
    --num-clusters 40 --num-restarts 10 --max-iters 10000 \
    --density binomial --seed 93521 > ${outdir}/pyclone_analysis/log.pyclone 
/mnt/ddngs/zhangsw/miniconda3/envs/pyclone-vi/bin/pyclone-vi write-results-file \
    -i ${outdir}/pyclone_analysis/${prefix}.pyclone_dat.h5 \
    -o ${outdir}/pyclone_analysis/${prefix}.pyclone_dat.h5.tsv

## old pyclone
#/mnt/ddngs/zhangsw/miniconda3/envs/pyclone/bin/PyClone run_analysis_pipeline \
#    --in_files ${outdir}/pyclone_analysis/${prefix}.pyclone_dat.tsv \
#    --working_dir ${outdir}/pyclone_analysis/ \
#    --burnin 1000 --num_iters 10000 --max_clusters 40 \
#    --density pyclone_binomial --seed 93521 > ${outdir}/pyclone_analysis/log.pyclone 2>&1

## site selection

python ${scriptDir}/site_selection_by_pyclone.py \
    --tsv ${outdir}/pyclone_analysis/${prefix}.pyclone_dat.h5.tsv \
    --dat ${outdir}/pyclone_analysis/${prefix}.pyclone_dat.tsv \
    --ref ${refFa} \
    --out ${outdir}/pyclone_analysis/${prefix}.pyclone_dat.site_selection.tsv \
    --hotspot /mnt/ddngs/zhangsw/project/MRD/OV_WES/hotspots/TCGA_merged.hotspot.tsv \
    --site 50 --gc 30,70 --flank 150 --min 3




## QC
#python2 mapping_QC_module_HR.siwen20211109.py \
#    --input_target ~{targetInterval}  \
#    --input_whole_region_info ~{coverage} \
#    --input_insert_size_info ~{insertsize} \
#    --output_file ~{outdir}/~{name}.qcStat.tsv \
#    --perdep ~{basedep} \
#    --bamfile ~{bam} \
#    --samtools ~{samtools}
