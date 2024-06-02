#!/usr/bin/bash

# <<<--- variables --->>>

shellname=$0
scriptDir=`dirname ${shellname}`
refFa='/mnt/user/errand/pipeline/MC3/ref/Homo_sapiens_assembly19.fasta'
#bed='/mnt/ddngs/zhangsw/project/MRD/OV_WES/bed/xgen-exome-research-panel-targets.b37.bed'
#interval='/mnt/user/errand/pipeline/MC3/ref/Agilent.SureSelect_All_Exon-v6_r2.interval_list'
annoDir='/mnt/ddngs/zhangsw/database/GRCh37'
snp_pileup='/mnt/user/errand/pipeline/OncoDeficiencyPro/tools/snppileup/htstools-master/snp-pileup'
#run_facet='/mnt/ddngs/lvr/script/HRD/OncoDeficiencyPro/script/run_facet.lr.R'
#facets_score='/mnt/ddngs/lvr/script/HRD/OncoDeficiencyPro/script/get_scarHRD_score_plus.R'
Rscript='/mnt/ddngs/lvr/Miniconda3/miniconda3/envs/py2.7/bin/Rscript'
#hotspot='/mnt/ddngs/zhangsw/project/MRD/OV_WES/hotspots/TCGA_merged.hotspot.tsv'
outdir=`pwd`'/outdir'
prefix='sample'
site=30

# usage
function usage {
    printf -- "Usage:\n    sh ${shellname} -b tumorBam -n normalBam -v somaticVcf -g germlineVcf -o outdir -p prefix -s 30 -l bed -t hotspot \n\n";
    printf -- "Required arguments: \n";
    printf -- "    -b|--tbam          tumor bam file \n";
    printf -- "    -n|--nbam          normal bam file \n";
    printf -- "    -v|--somaticVcf    tumor somatic vcf \n"
    printf -- "    -g|--germlineVcf   normal germline vcf \n";
    printf -- "Optional arguments: \n";
    printf -- "    -o|--outdir        output path, default '${outdir}' \n";
    printf -- "    -p|--prefix        output prefix, default '${prefix}' \n";
    printf -- "    -r|--ref           reference fasta, default '${refFa}' \n";
    printf -- "    -s|--site          site selection number, 30 for default \n";
    printf -- "    -l|--bed           bed region \n";
    printf -- "    -t|--hotspot       hotspot region, bed-like file, 
                         e.g /mnt/ddngs/zhangsw/project/MRD/OV_WES/hotspots/TCGA_merged.hotspot.tsv 
                         e.g /mnt/ddngs/zhangsw/project/MRD/OV_WES/hotspots/oncogene_whitelist.annotationCombined.filtered.tsv \n";
    printf -- "    -h|--help          help \n";
    printf -- "@Author: Siwen Zhang, zhang.siwen@puruijizhun.com \n\n";
    exit 0
}


# Arguments passing
TEMP=`getopt -o b:n:v:g:o:p:r:s:l:t:h --long tbam:,nbam:,somaticVcf:,germlineVcf:,outdir:,prefix:,ref:,site:,bed:,hotspot:,help \
     -n ${shellname} -- "$@"`

# no arguments detected
if [ $? != 0 ] || [ $# == 0 ]
    then usage
fi

# set arguments
eval set -- "$TEMP"
while true ; do
    case "$1" in
        -b|--tbam)
            tbam="$2"; shift 2;;
        -n|--nbam)
            nbam="$2"; shift 2;;
        -v|--somaticVcf)
            somaticVcf="$2"; shift 2;;
        -g|--germlineVcf)
            germlineVcf="$2"; shift 2 ;;
        -o|--outdir)
            outdir="$2"; shift 2 ;;
        -p|--prefix)
            prefix="$2"; shift 2 ;;
        -r|--ref)
            refFa="$2"; shift 2 ;;
        -s|--site)
            site="$2"; shift 2 ;;
        -l|--bed)
            bed="$2"; shift 2 ;;
        -t|--hotspot)
            hotspot="$2"; shift 2 ;;
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
if [ -f ${tbam} ] && [ -f ${nbam} ] && [ -f ${somaticVcf} ] && [ -f ${germlineVcf} ]
then echo "Files found."
    echo "Tumor BAM: ${tbam}"
    echo "Normal BAM: ${nbam}"
    echo "Somatic VCF: ${somaticVcf}"
    echo "Germline VCF: ${germlineVcf}"
else usage
fi



# <<<--- SNV indel calling --->>>

## filter blacklist (Duke, DAC, repeatMasker, simpleRepeat)

mkdir -p ${outdir}/snv_cnv/
if [ -f ${bed} ]
then
    bedtools subtract -nonamecheck -header -a ${somaticVcf} -b ${annoDir}/wgEncodeDacMapabilityConsensusExcludable.bed.gz |\
        bedtools subtract -nonamecheck -header -a - -b ${annoDir}/wgEncodeDukeMapabilityRegionsExcludable.bed.gz |\
        bedtools subtract -nonamecheck -header -a - -b ${annoDir}/genomicSuperDups.bed |\
        bedtools subtract -nonamecheck -header -a - -b ${annoDir}/hg19_repeatMasker.bed |\
        bedtools intersect -nonamecheck -header -u -a - -b ${bed} > ${outdir}/snv_cnv/${prefix}.somatic.filterBlacklist.vcf
    bcftools norm -O v -o ${outdir}/snv_cnv/${prefix}.somatic.filterBlacklist.norm.vcf -m - ${outdir}/snv_cnv/${prefix}.somatic.filterBlacklist.vcf
    #bcftools view --types snps ${outdir}/snv_cnv/${prefix}.somatic.filterBlacklist.norm.vcf > ${outdir}/snv_cnv/${prefix}.somatic.filterBlacklist.norm.onlysnp.vcf 

    bedtools subtract -nonamecheck -header -a ${germlineVcf} -b ${annoDir}/wgEncodeDacMapabilityConsensusExcludable.bed.gz |\
        bedtools subtract -nonamecheck -header -a - -b ${annoDir}/wgEncodeDukeMapabilityRegionsExcludable.bed.gz |\
        bedtools subtract -nonamecheck -header -a - -b ${annoDir}/genomicSuperDups.bed |\
        bedtools subtract -nonamecheck -header -a - -b ${annoDir}/hg19_repeatMasker.bed |\
        bedtools intersect -nonamecheck -header -u -a - -b ${bed} > ${outdir}/snv_cnv/${prefix}.germline.filterBlacklist.vcf
    bcftools norm -O v -o ${outdir}/snv_cnv/${prefix}.germline.filterBlacklist.norm.vcf -m - ${outdir}/snv_cnv/${prefix}.germline.filterBlacklist.vcf
else
    bedtools subtract -nonamecheck -header -a ${somaticVcf} -b ${annoDir}/wgEncodeDacMapabilityConsensusExcludable.bed.gz |\
        bedtools subtract -nonamecheck -header -a - -b ${annoDir}/wgEncodeDukeMapabilityRegionsExcludable.bed.gz |\
        bedtools subtract -nonamecheck -header -a - -b ${annoDir}/genomicSuperDups.bed |\
        bedtools subtract -nonamecheck -header -a - -b ${annoDir}/hg19_repeatMasker.bed > ${outdir}/snv_cnv/${prefix}.somatic.filterBlacklist.vcf
    bcftools norm -O v -o ${outdir}/snv_cnv/${prefix}.somatic.filterBlacklist.norm.vcf -m - ${outdir}/snv_cnv/${prefix}.somatic.filterBlacklist.vcf
    #bcftools view --types snps ${outdir}/snv_cnv/${prefix}.somatic.filterBlacklist.norm.vcf > ${outdir}/snv_cnv/${prefix}.somatic.filterBlacklist.norm.onlysnp.vcf 

    bedtools subtract -nonamecheck -header -a ${germlineVcf} -b ${annoDir}/wgEncodeDacMapabilityConsensusExcludable.bed.gz |\
        bedtools subtract -nonamecheck -header -a - -b ${annoDir}/wgEncodeDukeMapabilityRegionsExcludable.bed.gz |\
        bedtools subtract -nonamecheck -header -a - -b ${annoDir}/genomicSuperDups.bed |\
        bedtools subtract -nonamecheck -header -a - -b ${annoDir}/hg19_repeatMasker.bed > ${outdir}/snv_cnv/${prefix}.germline.filterBlacklist.vcf
    bcftools norm -O v -o ${outdir}/snv_cnv/${prefix}.germline.filterBlacklist.norm.vcf -m - ${outdir}/snv_cnv/${prefix}.germline.filterBlacklist.vcf
fi



# <<<--- Copy Number Variation --->>>

## run facets

mkdir -p ${outdir}/facets
${snp_pileup} -g -q1 -Q20 ${outdir}/snv_cnv/${prefix}.germline.filterBlacklist.norm.vcf ${outdir}/facets/${prefix}.snp-pileup ${nbam} ${tbam}

${Rscript} ${scriptDir}/run_facet.lr.zsw.R ${outdir}/facets/${prefix}.snp-pileup.gz ${prefix} ${outdir}/facets/${prefix}



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

python ${scriptDir}/pyclone_dat_generation.v2.py \
    --infile ${outdir}/snv_cnv/${prefix}.somatic.filterBlacklist.norm.maf \
    --cnv ${outdir}/facets/${prefix}/${prefix}.4.segment.final.txt \
    --out ${outdir}/pyclone_analysis/${prefix}.pyclone_dat.tsv \
    --gender female --sampleid ${prefix} \
    --depth 0 --ad 0 --af 0 --source facets --indel 5 

## run new pyclone

/mnt/ddngs/zhangsw/miniconda3/envs/pyclone-vi/bin/pyclone-vi fit \
    -i ${outdir}/pyclone_analysis/${prefix}.pyclone_dat.tsv \
    -o ${outdir}/pyclone_analysis/${prefix}.pyclone_dat.h5 \
    --num-clusters 40 --num-restarts 10 --max-iters 10000 \
    --density binomial --seed 93521 > ${outdir}/pyclone_analysis/log.pyclone 
/mnt/ddngs/zhangsw/miniconda3/envs/pyclone-vi/bin/pyclone-vi write-results-file \
    -i ${outdir}/pyclone_analysis/${prefix}.pyclone_dat.h5 \
    -o ${outdir}/pyclone_analysis/${prefix}.pyclone_result.tsv


## site selection

if [ "${hotspot}" == ""]
then
    python ${scriptDir}/site_selection_by_pyclone.v2.py \
        --tsv ${outdir}/pyclone_analysis/${prefix}.pyclone_result.tsv \
        --dat ${outdir}/pyclone_analysis/${prefix}.pyclone_dat.tsv \
        --ref ${refFa} \
        --out ${outdir}/pyclone_analysis/${prefix}.site_selection${site}.tsv \
        --maf ${outdir}/snv_cnv/${prefix}.somatic.filterBlacklist.norm.maf \
        --site ${site} --gc 30,70 --flank 100 --min 3
else
    python ${scriptDir}/site_selection_by_pyclone.v2.py \
        --tsv ${outdir}/pyclone_analysis/${prefix}.pyclone_result.tsv \
        --dat ${outdir}/pyclone_analysis/${prefix}.pyclone_dat.tsv \
        --ref ${refFa} \
        --out ${outdir}/pyclone_analysis/${prefix}.site_selection${site}.tsv \
        --hotspot ${hotspot} \
        --maf ${outdir}/snv_cnv/${prefix}.somatic.filterBlacklist.norm.maf \
        --site ${site} --gc 30,70 --flank 100 --min 3
fi



## QC
#python2 mapping_QC_module_HR.siwen20211109.py \
#    --input_target ~{targetInterval}  \
#    --input_whole_region_info ~{coverage} \
#    --input_insert_size_info ~{insertsize} \
#    --output_file ~{outdir}/~{name}.qcStat.tsv \
#    --perdep ~{basedep} \
#    --bamfile ~{bam} \
#    --samtools ~{samtools}
