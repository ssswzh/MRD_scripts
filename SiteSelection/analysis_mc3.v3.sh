#!/usr/bin/bash

# <<<--- variables --->>>

shellname=$0
scriptDir=`dirname ${shellname}`
refFa='/mnt/user/errand/pipeline/MC3/ref/Homo_sapiens_assembly19.fasta'
hotspot='/mnt/ddngs/zhangsw/project/MRD/hotspots/oncogene_whitelist.annotationCombined.filtered.tsv'
#bed='/mnt/ddngs/zhangsw/project/MRD/OV_WES/bed/xgen-exome-research-panel-targets.b37.bed'
#interval='/mnt/user/errand/pipeline/MC3/ref/Agilent.SureSelect_All_Exon-v6_r2.interval_list'
pycloneVi=/mnt/ddngs/zhangsw/miniconda3/envs/pyclone-vi/bin/pyclone-vi
annoDir='/mnt/ddngs/zhangsw/database/GRCh37'
snp_pileup='/mnt/user/errand/pipeline/OncoDeficiencyPro/tools/snppileup/htstools-master/snp-pileup'
Rscript='/mnt/ddngs/lvr/Miniconda3/miniconda3/envs/py2.7/bin/Rscript'
python=/mnt/ddngs/zhangsw/miniconda3/bin/python
outdir=`pwd`'/outdir'
prefix='sample'
site=25
bed=/mnt/ddngs/zhangsw/project/MRD/bed/xgen-exome-research-panel-targets.b37.bed
purity=""

# usage
function usage {
    printf -- "Usage:\n    sh ${shellname} -i /path/to/OncoBuster/result/ -o outdir -p prefix -s 25 -l bed -t hotspot.tsv \n\n";
    printf -- "Required arguments: \n";
    printf -- "    -i|--indir         output directory of OncoBuster 
                        e.g. /mnt/ddngs/errand/Project/OncoBuster/Clinical/LJ2205318/ \n";
    printf -- "Optional arguments: \n";
    printf -- "    -o|--outdir        output path, default '${outdir}' \n";
    printf -- "    -p|--prefix        output prefix, default '${prefix}' \n";
    printf -- "    -r|--ref           reference fasta, default '${refFa}' \n";
    printf -- "    -s|--site          site selection number, default '${site}' \n";
    printf -- "    -l|--bed           bed region, default '${bed}' \n";
    printf -- "    -t|--hotspot       hotspot region, bed-like file, default '${hotspot}'
                         e.g /mnt/ddngs/zhangsw/project/MRD/hotspots/TCGA_merged.hotspot.tsv 
                         e.g /mnt/ddngs/zhangsw/project/MRD/hotspots/oncogene_whitelist.annotationCombined.filtered.tsv \n";
    printf -- "    -y|--purity        purity from clinical image estimation, default NULL \n";
    printf -- "    -h|--help          help \n";
    printf -- "@Author: Siwen Zhang, zhang.siwen@puruijizhun.com \n\n";
    exit 1
}


# Arguments passing
TEMP=`getopt -o i:o:p:r:s:l:t:y:h --long indir:,outdir:,prefix:,ref:,site:,bed:,hotspot:,purity:,help \
     -n ${shellname} -- "$@"`

# no arguments detected
if [ $? != 0 ] || [ $# == 0 ]
    then usage
fi

# set arguments
eval set -- "$TEMP"
while true ; do
    case "$1" in
        -i|--indir)
            indir="$2"; shift 2;;
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
        -y|--purity)
            purity="$2"; shift 2 ;;
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

# files not found
tbam=${indir}/Mapping/Tumor.final.bam
nbam=${indir}/Mapping/Normal.final.bam
somaticVcf=${indir}/SomaticMutation/MC3_snake/merge_vcfs.roi.vep.vcf
germlineVcf=${indir}/GermlineMutation/germline.final.vcf
if [ -f ${tbam} ] && [ -f ${nbam} ] && [ -f ${somaticVcf} ] && [ -f ${germlineVcf} ]
then echo "Files found."
    echo "Tumor BAM: ${tbam}"
    echo "Normal BAM: ${nbam}"
    echo "Somatic VCF: ${somaticVcf}"
    echo "Germline VCF: ${germlineVcf}"
else echo "Files not found: " 
    echo "Tumor BAM: ${tbam}"
    echo "Normal BAM: ${nbam}"
    echo "Somatic VCF: ${somaticVcf}"
    echo "Germline VCF: ${germlineVcf}"
    usage
fi



# <<<--- SNV indel calling --->>>

## filter blacklist (Duke, DAC, repeatMasker, simpleRepeat)

# filter somatic vcf, at least two softwares dectected the same site

mkdir -p ${outdir}/snv_cnv/

${python} ${scriptDir}/src/filter_caller.py --vcf ${somaticVcf} --out ${outdir}/snv_cnv/merge_vcfs.filter_softwares.vcf
somaticVcf=${outdir}/snv_cnv/merge_vcfs.filter_softwares.vcf


if [  "${bed}" != "" ]
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

${Rscript} ${scriptDir}/src/run_facet.lr.zsw.R ${outdir}/facets/${prefix}.snp-pileup.gz ${prefix} ${outdir}/facets/${prefix}



# <<<--- PyClone cluster and phylogenetic tree construction --->>>

## prepare data

mkdir -p ${outdir}/pyclone_analysis/

perl ${scriptDir}/src/vcf2maf.pl \
    --input-vcf ${outdir}/snv_cnv/${prefix}.somatic.filterBlacklist.norm.vcf \
    --output-maf ${outdir}/snv_cnv/${prefix}.somatic.filterBlacklist.norm.maf \
    --tumor-id PRIMARY --normal-id NORMAL \
    --custom-enst ${scriptDir}/db/GCF_000001405.25_GRCh37.p13_genomic.selected_transcript.txt \
    --inhibit-vep --ref-fasta ${refFa} \
    --filter-vcf 0 --ncbi-build GRCh37 --retain-info COSMIC,CENTERS,CONTEXT,DBVS


if [ "${purity}" == "" ]
then known_purity=""
else known_purity="--purity ${purity}"
fi

${python} ${scriptDir}/src/pyclone_dat_generation.v2.py \
    --infile ${outdir}/snv_cnv/${prefix}.somatic.filterBlacklist.norm.maf \
    --cnv ${outdir}/facets/${prefix}/${prefix}.4.segment.final.txt \
    --out ${outdir}/pyclone_analysis/${prefix}.pyclone_dat.tsv \
    --gender female --sampleid ${prefix} \
    --depth 0 --ad 0 --af 0 --source facets ${known_purity}


## run new pyclone

${pycloneVi} fit \
    -i ${outdir}/pyclone_analysis/${prefix}.pyclone_dat.tsv \
    -o ${outdir}/pyclone_analysis/${prefix}.pyclone_dat.h5 \
    --num-clusters 40 --num-restarts 10 --max-iters 10000 \
    --density binomial --seed 93521 > ${outdir}/pyclone_analysis/log.pyclone 
${pycloneVi} write-results-file \
    -i ${outdir}/pyclone_analysis/${prefix}.pyclone_dat.h5 \
    -o ${outdir}/pyclone_analysis/${prefix}.pyclone_result.tsv


## site selection

if [ "${hotspot}" == "" ]
then known_hotspot=""
else known_hotspot="--hotspot ${hotspot} "
fi

${python} ${scriptDir}/src/site_selection_by_pyclone.v3.py \
    --tsv ${outdir}/pyclone_analysis/${prefix}.pyclone_result.tsv \
    --dat ${outdir}/pyclone_analysis/${prefix}.pyclone_dat.tsv \
    --ref ${refFa} \
    --out ${outdir}/pyclone_analysis/${prefix}.site_selection${site}.tsv \
    --maf ${outdir}/snv_cnv/${prefix}.somatic.filterBlacklist.norm.maf \
    --site ${site} --af 0.05 --gc 30,70 --flank 100 --min 3 --prev AF:0.01 --indel 8,30 --func high ${known_hotspot}

cut -f21,25-27,55-56 ${outdir}/pyclone_analysis/${prefix}.site_selection${site}.tsv|sed '1d' > ${outdir}/pyclone_analysis/${prefix}.site_selection${site}.order.tsv

${python} ${scriptDir}/src/fill_order_template.py --tsv ${outdir}/pyclone_analysis/${prefix}.site_selection${site}.order.tsv --id ${prefix} --out ${outdir}/${prefix}.site_selection${site}.order.xlsx --temp ${scriptDir}/db/MRDTargetDesign.PillarBiosciences.template.xlsx
