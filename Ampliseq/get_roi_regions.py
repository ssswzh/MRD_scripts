#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Author  : zhangsiwen
# @ChangeLog:
#       20230214, first version
#       20230315, change site output to bed format, start = start -1
#       20230426, add MultiRegionSplitter() to check multi site in single region and expand sites


import argparse
import pandas as pd


def GetArgs():
    parser = argparse.ArgumentParser(description='extract amplicon ROI region and site from design', formatter_class=argparse.RawTextHelpFormatter)
    parser._action_groups.pop()
    # required arguments
    required = parser.add_argument_group('Required arguments')
    required.add_argument('--roi', help='amplicon ROI region design by Zhengu in XLSX format, \nheaders are: Target_Name\tTarget_Coordinates\tChrom\tROC_Start\tROC_End\tGene_Symbol', action='store', dest='roi', required=True)
    required.add_argument('--site', help='site selection result from MRD SiteSelection final result, endswith site_selection{number}.tsv', action='store', dest='site', required=True)
    required.add_argument('--out', help='out file prefix, two files will be generated, endswith site_selection.final.bed and roi_region.bed separately', action='store', dest='out', required=True)
    usage = '''Usage:
    %(prog)s --design purui_{sample}_external_primer_design.xlsx --site {sample}.site_selection25.tsv --out outprefix
'''
    parser.epilog = usage
    args = parser.parse_args()
    return args


'''
design = 'purui_990370002_external_primer_design_20221229.xlsx'
site = '990370002.site_selection25.tsv'
out = 'test'
'''


def getHeader(column_names):
    if 'Chrom' in column_names:
        return 'Chrom', 'ROC_Start', 'ROC_End', 'Target_Name'
    if '染色体' in column_names:
        return '染色体', 'ROI起始位置', 'ROI终止位置', '名称'


def MultiRegionSplitter(df, col):
    column_order = df.columns
    df = df.drop(col, axis=1) \
            .join(df[col] \
            .str.split('\|\|', expand=True) \
            .stack() \
            .reset_index(level=1, drop=True).rename(col)) \
            .reset_index(drop=True)
    df = df[df.columns]
    return df


def main():
    ''' '''
    args = GetArgs()
    roi_design = pd.read_excel(args.roi)
    site_design = pd.read_csv(args.site, sep='\t')
    ROCchrom, ROCstart, ROCend, targetName = getHeader(roi_design.columns)
    roi_design[ROCchrom] = roi_design[ROCchrom].str.replace('chr','')
    roi_design[ROCstart] = roi_design[ROCstart].astype(int)
    roi_design[ROCend] = roi_design[ROCend].astype(int)

    if [i for i in roi_design[targetName] if "||" in i]:
        roi_design = MultiRegionSplitter(df=roi_design, col=targetName)

    roi_design[['index', 'Gene']] = roi_design[targetName].str.split('.', 1, expand=True)
    roi_design['index'] = roi_design['index'].str.replace('ZG','').astype(int)
    roi_design = roi_design.sort_values('index')
    site_design['index'] = [i+1 for i in site_design.index.to_list()]

    combinedDf = pd.merge(roi_design, site_design, left_on='index', right_on='index', how='inner')
    if 'Gene' in combinedDf.columns:
        geneCol = 'Gene'
    else:
        geneCol = 'Gene_x'
    roi_intersect = combinedDf[[ROCchrom, ROCstart, ROCend,'index',geneCol,'mutation_id']]
    roi_intersect.to_csv(args.out+".roi_region.final.bed", sep='\t', index=False, mode='w', header=False)

    site_cols = ['Chromosome','Start_Position','End_Position','Reference_Allele','Tumor_Seq_Allele2','vaf','Hugo_Symbol','Variant_Classification']
    site_intersect = combinedDf[site_cols]
    site_intersect.columns = str('#' + ','.join(site_cols)).split(',')
    site_intersect['Start_Position'] = site_intersect['Start_Position'] - 1
    site_intersect.to_csv(args.out+".site_selection.final.bed", sep='\t', index=False, mode='w', header=True)

    return


if __name__ == '__main__':
    main()
