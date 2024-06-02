# MRD分析流程

MRD项目的分析流程，包括选点分析和扩增子分析两部分。

## 扩增子变异检测 Variants Detection From Amplicon Sequencing Data

### 检测原理

MRD扩增子分析要求对一个患者的样本（检测样本）配置至少5个阴性对照样本（此处为2022开头的健康人和NA开头的GIAB标准品），患者样本可以接收多个，例如[示例文件](###_流程使用示例_script_usage_example)中提供了相同患者的两个检测样本。

对照样本的作用是对每个位点进行[polishing](https://git.puruijizhun.com/siwen.zhang/mrdanalysis/ampliseq/src/polishing_rebuild.r)，过滤掉测序错误导致的超低频突变。

首先对每个样本（不区分对照样本和检测样本）进行突变检测，等到一个批次的样本都检测完后，再合并一个批次样本的深度、频率，然后进行位点抛光，得到最终检测样本的突变阴阳性。


### 软件配置 Software Configuration

流程需要的软件有以下四个，均写在单样本分析流程`analysis_pipeline.v3.sh`中：

```
VardictBin=/mnt/ddngs/zhangsw/software/VarDict-1.8.3/bin
fastp=/mnt/ddngs/zhangsw/bin/fastp
picard=/mnt/ddngs/zhangsw/miniconda3/bin/picard
Rscript='/mnt/ddngs/zhangsw/miniconda3/bin/Rscript'
```

其他常用软件未进行配置，均可在`/mnt/user/share/bin/`中找到，包括`usearch`，`sambamba`，`bwa`，`bcftools`等。

### 流程使用说明 Script Usage Instruction

本分析流程采用bash、python、R语言进行编写，流程脚本储存位置为`/mnt/ddngs/zhangsw/project/MRD/scripts/Ampliseq/`，共包含1个多样本批量主分析脚本、1个单样本分析脚本、7个不同功能的脚本，总体分析脚本集成各功能模块于`analysis_ampliseq.batch.v1.sh`，使用方法如下：

```
Usage:
    sh analysis_ampliseq.batch.v1.sh -i config -l bed -s site.bed -v refFa -o Outdir -p prefix -t thread -c low

Required arguments:
    -i|--input         input config file, tab-delimited, Group column only accept 'control' or 'test',
                       Sample column will be used as output dir name and file prefix name
                       Sample name CAN NOT start with number!!!
    -l|--bed           bed region
    -s|--site          site coordinate bed
Optional arguments:
    -v|--refFa         reference fasta, default '/mnt/ddngs/zhangsw/database/GRCh37_clinical/nochrY/Homo_sapiens_assembly19.nochrY.fasta'
    -o|--outdir        output path, default '/mnt/ddngs/zhangsw/project/MRD/scripts/Ampliseq/Outdir'
    -t|--thread        thread, default ''
    -c|--clean         if clean intermediate files, levels: none|low|high, default 'low'
                       none: keep all results,
                       low: remove fastp/*fq, fastq_merge/*merged.fq, mapping/*sam, mapping/unsorted.bam,
                       high: remove fastp/, fastq_merge/, mapping/,
    -h|--help          help

Input file example, sample name can not start with number:
    Sample  R1  R2  Group
    N1  N1_R1.fq    N1_R2.fq    control
    N2  N2_R1.fq    N2_R2.fq    control
    N3  N3_R2.fq    N3_R3.fq    control
    N4  N4_R2.fq    N4_R3.fq    control
    N5  N5_R3.fq    N5_R4.fq    control
    S1  S1_R1.fq    S1_R2.fq    test
    S2  S2_R1.fq    S2_R2.fq    test

@Author: Siwen Zhang

```

另外，对单样本检测突变的脚本`analysis_pipeline.v3.sh`会被主分析脚本`analysis_ampliseq.batch.v1.sh`调用。

### 流程使用示例 Script Usage Example

流程脚本储存位置为`${scriptDir}`，测试目录为`/mnt/ddngs/zhangsw/project/MRD/test_pipeline/test_Ampliseq/`，`sample_list`内容如下，**要求第一列不能以数字开头**，有无表头均可：

| Sample | R1       | R2       | Group   |
|--------|----------|----------|---------|
| N20220623005BC03-6 | /mnt/ddngs/zhangsw/project/MRD/test_pipeline/test_Ampliseq/data/20220623005BC03-6_1.fq.gz | /mnt/ddngs/zhangsw/project/MRD/test_pipeline/test_Ampliseq/data/20220623005BC03-6_2.fq.gz | control |
| S990370003BC01-8   | /mnt/ddngs/zhangsw/project/MRD/test_pipeline/test_Ampliseq/data/990370003BC01-8_1.fq.gz   | /mnt/ddngs/zhangsw/project/MRD/test_pipeline/test_Ampliseq/data/990370003BC01-8_2.fq.gz   | test    |
| LJ2103907BC01-7   | /mnt/ddngs/zhangsw/project/MRD/test_pipeline/test_Ampliseq/data/LJ2103907BC01-7_1.fq.gz   | /mnt/ddngs/zhangsw/project/MRD/test_pipeline/test_Ampliseq/data/LJ2103907BC01-7_2.fq.gz   | test    |
| NA12878-1         | /mnt/ddngs/zhangsw/project/MRD/test_pipeline/test_Ampliseq/data/NA12878-1_1.fq.gz         | /mnt/ddngs/zhangsw/project/MRD/test_pipeline/test_Ampliseq/data/NA12878-1_2.fq.gz         | control |
| NA24143-3         | /mnt/ddngs/zhangsw/project/MRD/test_pipeline/test_Ampliseq/data/NA24143-3_1.fq.gz         | /mnt/ddngs/zhangsw/project/MRD/test_pipeline/test_Ampliseq/data/NA24143-3_2.fq.gz         | control |
| NA24149-2         | /mnt/ddngs/zhangsw/project/MRD/test_pipeline/test_Ampliseq/data/NA24149-2_1.fq.gz         | /mnt/ddngs/zhangsw/project/MRD/test_pipeline/test_Ampliseq/data/NA24149-2_2.fq.gz         | control |
| NA24694-4         | /mnt/ddngs/zhangsw/project/MRD/test_pipeline/test_Ampliseq/data/NA24694-4_1.fq.gz         | /mnt/ddngs/zhangsw/project/MRD/test_pipeline/test_Ampliseq/data/NA24694-4_2.fq.gz         | control |
| NA24695-5         | /mnt/ddngs/zhangsw/project/MRD/test_pipeline/test_Ampliseq/data/NA24695-5_1.fq.gz         | /mnt/ddngs/zhangsw/project/MRD/test_pipeline/test_Ampliseq/data/NA24695-5_2.fq.gz         | control |

**注意**bed文件的区间是不包含扩增引物的区间，即`ROI_Start`和`ROI_End`，不是`Amplicon_Start`和`Amplicon_End`。

分析代码示例如下：

```
# 启动分析
source /mnt/ddngs/zhangsw/.bashrc
sh ${scriptDir}/analysis_ampliseq.batch.v1.sh \
    -i sample_list \
    -l bed/SAM001_internal_primer_design_20220822.bed 
    -s bed/SAM001.site_selection.bed \
    -v /mnt/ddngs/zhangsw/database/GRCh37_clinical/nochrY/Homo_sapiens_assembly19.nochrY.fasta \
    -o analysis \
    -t 8 -c low
```

**可将上述代码写入`run.sh`，然后直接运行`sh run.sh > run.sh.log`，`analysis_ampliseq.batch.v1.sh`会自动使用`qsub`投递任务。**


***

## 扩增子变异检测最终结果 Final Variant Results From Amplicon Sequencing Data

分析流程主要的结果有两类，一类是每个样本的数据质控，一类是检测样本的变异阴阳性判定。

### 样本数据质控 Sequencing Data Quality Control

数据质控有两个主要结果，一个是样本整体比对水平，一个是每个区间深度分布。

**样本整体比对水平**

`input`文件中的各个样本的质控结果储存在`${outdir}/[control/test]-{sample}/{sample}.all.QCstats.tsv`中。

**每个区间深度分布**

每个区间的深度分布呈现格式为图片PDF，储存在`${outdir}/polishing/allsample.DP.boxplot.pdf`，是所有样本的区间深度的箱线图。

### 检测样本的变异阴阳性判定 Determination of The Variant Positive of The Test Sample

结果储存在`${outdir}/polishing/all_sample.final.result.tsv`中，为TAB分割的文本文件，前21列为阴性对照样本对位点的检测结果，这里不进行展示。从第22列开始是检测样本的位点检测结果，具体结果如下表：

| chrom | pos       | ref | alt | LJ2103907BC01.7 | S990370003BC01.8 | NoiseType | id          | Result:LJ2103907BC01.7 | Result:S990370003BC01.8 |
|-------|-----------|-----|-----|-----------------|------------------|-----------|-------------|------------------------|-------------------------|
| 15    | 78393788  | C   | T   | 0.000116219     | 6.61367e-05      | C>T       | 15:78393788 | Positive               | Negative                |
| 14    | 91809981  | T   | A   | 0               | 0                | T>A       | 14:91809981 | Negative               | Negative                |
| 6     | 157469867 | A   | G   | 0.000159648     | 9.94588e-05      | A>G       | 6:157469867 | Positive               | Negative                |
| 6     | 159653030 | T   | C   | 0.000318086     | 0.000152478      | T>C       | 6:159653030 | Positive               | Negative                |
| 11    | 6643140   | G   | T   | 0               | 0                | G>T       | 11:6643140  | Negative               | Negative                |
| 9     | 135940607 | G   | A   | 0.000327938     | 1.53947e-05      | G>A       | 9:135940607 | Positive               | Negative                |
| 1     | 5965418   | G   | C   | 9.17291e-05     | 0                | G>C       | 1:5965418   | Positive               | Negative                |
| 1     | 159558339 | A   | G   | 0.000361071     | 0.000178859      | A>G       | 1:159558339 | Positive               | Negative                |
| 22    | 19966530  | G   | C   | 0.000146759     | 0                | G>C       | 22:19966530 | Positive               | Negative                |
| 10    | 82348467  | A   | T   | 0               | 0                | A>T       | 10:82348467 | Negative               | Negative                |
| 10    | 26800679  | G   | A   | 7.24423e-05     | 7.95097e-05      | G>A       | 10:26800679 | Negative               | Negative                |
| 2     | 233407227 | A   | T   | 0               | 0                | A>T       | 2:233407227 | Negative               | Negative                |
| 19    | 2137738   | C   | G   | 0.000193792     | 0                | C>G       | 19:2137738  | Positive               | Negative                |
| 1     | 236205571 | C   | A   | 0               | 0                | C>A       | 1:236205571 | Negative               | Negative                |
| 15    | 42167140  | T   | A   | 0               | 0                | T>A       | 15:42167140 | Negative               | Negative                |
| 10    | 14563963  | G   | A   | 0.000223146     | 0.000163203      | G>A       | 10:14563963 | Negative               | Negative                |
| X     | 153229648 | G   | A   | 0               | 2.05808e-05      | G>A       | X:153229648 | Negative               | Negative                |
| 6     | 35103933  | C   | G   | 0.000244606     | 0                | C>G       | 6:35103933  | Positive               | Negative                |
| 11    | 94731606  | A   | T   | 0.000140508     | 0                | A>T       | 11:94731606 | Positive               | Negative                |
| 12    | 81719624  | T   | A   | 0.000153192     | 0                | T>A       | 12:81719624 | Positive               | Negative                |

其中样本名称`{sample}`列为该样本实测的突变频率，`Result:{sample}`列为该样本在每个位点处判定的阴阳性结果。

**样本水平阴阳性判定规则为2个位点检出阳性。**
