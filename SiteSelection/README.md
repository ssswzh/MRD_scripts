# MRD分析流程

MRD项目的分析流程，包括选点分析和扩增子分析两部分。

## WES选点分析 Site Selection Analysis From Whole Exome Sequencing Data

### 软件配置 Software Configuration

流程需要的软件有以下四个，均写在主分析流程`analysis_mc3.v3.sh`中：

```
pycloneVi=/mnt/ddngs/zhangsw/miniconda3/envs/pyclone-vi/bin/pyclone-vi
annoDir='/mnt/ddngs/zhangsw/database/GRCh37'
snp_pileup='/mnt/user/errand/pipeline/OncoDeficiencyPro/tools/snppileup/htstools-master/snp-pileup'
Rscript='/mnt/ddngs/lvr/Miniconda3/miniconda3/envs/py2.7/bin/Rscript'
```

其中，`annoDir`包含重复区间和不可比对区域：

- wgEncodeDacMapabilityConsensusExcludable.bed.gz
- wgEncodeDukeMapabilityRegionsExcludable.bed.gz
- genomicSuperDups.bed
- hg19_repeatMasker.bed

### 流程使用说明 Script Usage Instruction

本分析流程采用bash、python、R语言进行编写，流程脚本储存位置为`/mnt/ddngs/zhangsw/project/MRD/scripts/SiteSelection/`，共包含1个主分析脚本和11个不同功能的脚本，总体分析脚本集成各功能模块于`analysis_mc3.v3.sh`，使用方法如下：

```
Usage:
    sh analysis_mc3.v3.sh -i /path/to/OncoBuster/result/ -o outdir -p prefix -s 25 -l bed -t hotspot.tsv

Required arguments:
    -i|--indir         output directory of OncoBuster
                        e.g. /mnt/ddngs/errand/Project/OncoBuster/Clinical/LJ2205318/
Optional arguments:
    -o|--outdir        output path, default '/mnt/ddngs/zhangsw/project/MRD/scripts/SiteSelection/outdir'
    -p|--prefix        output prefix, default 'sample'
    -r|--ref           reference fasta, default '/mnt/user/errand/pipeline/MC3/ref/Homo_sapiens_assembly19.fasta'
    -s|--site          site selection number, default '25'
    -l|--bed           bed region, default '/mnt/ddngs/zhangsw/project/MRD/bed/xgen-exome-research-panel-targets.b37.flank20bp.bed'
    -t|--hotspot       hotspot region, bed-like file, default '/mnt/ddngs/zhangsw/project/MRD/OV_WES/hotspots/oncogene_whitelist.annotationCombined.filtered.tsv'
                         e.g /mnt/ddngs/zhangsw/project/MRD/OV_WES/hotspots/TCGA_merged.hotspot.tsv
                         e.g /mnt/ddngs/zhangsw/project/MRD/OV_WES/hotspots/oncogene_whitelist.annotationCombined.filtered.tsv
    -y|--purity        purity from clinical image estimation, default NULL
    -h|--help          help
@Author: Siwen Zhang, zhang.siwen@puruijizhun.com

```

### 流程使用示例 Script Usage Example

- 流程脚本储存位置为`${scriptDir}`
- 分析所需昂可星结果路径为`/mnt/ddngs/errand/Project/OncoBuster/Clinical/LJ2205318/`
- 分析所需参考基因组记为`/mnt/user/errand/pipeline/MC3/ref/Homo_sapiens_assembly19.fasta`
- 分析结果输出目录记为`/mnt/ddngs/zhangsw/project/MRD/test_pipeline/test_SiteSelection`
- 分析结果输出前缀记为`LJ2205318`
- 选点数量为`25`
- 包含的变异区间为昂可星flank20bp区间`/mnt/ddngs/zhangsw/project/MRD/bed/xgen-exome-research-panel-targets.b37.flank20bp.bed`
- 热点突变为`/mnt/ddngs/zhangsw/project/MRD/hotspots/oncogene_whitelist.annotationCombined.filtered.tsv`

分析代码示例如下：

```
# 启动分析
source /mnt/ddngs/zhangsw/.bashrc
sh ${scriptDir/analysis_mc3.v3.sh \
    --indir /mnt/ddngs/errand/Project/OncoBuster/Clinical/LJ2205318/ \
    --outdir /mnt/ddngs/zhangsw/project/MRD/test_pipeline/test_SiteSelection \
    --prefix LJ2205318 \
    --ref /mnt/user/errand/pipeline/MC3/ref/Homo_sapiens_assembly19.fasta \
    --site 25 \
    --bed /mnt/ddngs/zhangsw/project/MRD/bed/xgen-exome-research-panel-targets.b37.flank20bp.bed \
    --hotspot /mnt/ddngs/zhangsw/project/MRD/hotspots/oncogene_whitelist.annotationCombined.filtered.tsv
```

可将上述代码写入`run.sh`，然后使用`qsub`投递任务：

```
qsub -cwd -l mf=12g -q clinical.q@cn-thin05 run.sh
```

示例结果可见`/mnt/ddngs/zhangsw/project/MRD/test_pipeline/test_SiteSelection`，选点结果文件为`${outdir}/${sample}.site_selection${site}.xlsx`
