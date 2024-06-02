#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2022/01/13
# @Author  : zhangsiwen
# @contact : zhang.siwen@puruijizhun.com
# @LastModified: 2022/01/19
# @ChangeLog
#     20220113, first version
#     20220114, add --sort and --mode, new method for sorting sites
#     20220119, add Covered_Sample and Covered_Sample_Number column (sums) for new merged dataframe


import pandas as pd
import argparse
import re


def GetArgs():
    parser = argparse.ArgumentParser(description='Merge dataframes by columns (like cbind).', formatter_class=argparse.RawTextHelpFormatter)
    parser._action_groups.pop()
    required = parser.add_argument_group('Required arguments')
    required.add_argument('--files', help='input files, sepratated by comma ","', action='store', dest='files', required=True)
    required.add_argument('--out', help='output file name', action='store', dest='out', required=True)
    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument('--sort', help='sort by condition, default [CoveredSampleSum]\nCoveredSampleSum: sum of covered sample numbers across datasets\nAverageCoveredSamplePercent: average of covered sample percent across datasets', choices=['CoveredSampleSum','AverageCoveredSamplePercent'], action='store', default='CoveredSampleSum', dest='sort', required=False)
    optional.add_argument('--mode', help='merge mode, default "union"', choices=['intersect','union'], action='store', default='union', dest='mode', required=False)
    usage = '''Other info:
    Merge by columns:
        ['Chromosome','Start_Position','End_Position','Variant_Classification','Variant_Type','Reference_Allele']
    Rename columns by pattern:
        ['Occurrence','Covered_Sample', 'Covered_Sample_Number', 'Cumulate_Sample_Number', 'Cumulate_Sample_Number_Percent']
    '''
    parser.epilog = usage
    args = parser.parse_args()
    return args
    
'''
files = 'BRCA/brca_ink4_msk_2021/analysis/data_mutation.Summary.BRCA.HotSites.tsv,BRCA/brca_metabric/analysis/data_mutation.Summary.BRCA.HotSites.tsv,BRCA/brca_msk_2018/analysis/data_mutation.Summary.BRCA.HotSites.tsv,BRCA/brca_tcga/analysis/data_mutation.Summary.BRCA.HotSites.tsv,BRCA/brca_tcga_pan_can_atlas_2018/analysis/data_mutation.Summary.BRCA.HotSites.tsv,BRCA/brca_tcga_pub2015/analysis/data_mutation.Summary.BRCA.HotSites.tsv,BRCA/brca_tcga_pub/analysis/data_mutation.Summary.BRCA.HotSites.tsv'
files = 'TCGA.COAD.HotSites.tsv,TCGA.READ.HotSites.tsv'
sort = 'CoveredSampleSum'
mode = 'outer'
'''
cols = ['Chromosome','Start_Position','End_Position','Variant_Classification','Variant_Type']
pattern = ['Occurrence','Covered_Sample', 'Covered_Sample_Number', 'Covered_Sample_Number_Percent', 'Cumulate_Sample_Number', 'Cumulate_Sample_Number_Percent']


def RenameDF(df, prefix):
    new_column = []
    for i in list(df.columns):
        if i.startswith(tuple(pattern)):
            new_column.append(prefix+':'+i)
        else:
            new_column.append(i)
    return new_column


def AddColumnPercent(df, total, col, name):
    '''
    @ func: add percent column to dataframe
    '''
    # calculate cumulate sample percent for each site
    df[name] = df[col]/total
    return df


def ColumnsByPattern(columns, pattern):
    '''
    @ func: give a list, return elements with pattern
    '''
    return [i for i in list(columns) if i.endswith(pattern)]
    

def MergeDataframes(files, out, mode, sort):
    '''
    @ func: merge dataframes and sort
    '''
    files = [i.strip() for i in files.strip().split(',')]
    # read and merge dataframes
    for i in range(0,len(files)):
        if i == 0:
            mergedDf = pd.read_table(files[i], low_memory=False, comment='#')
            mergedDf = AddColumnPercent(mergedDf, max(mergedDf['Cumulate_Sample_Number']), 'Covered_Sample_Number', 'Covered_Sample_Number_Percent')
            mergedDf.columns = RenameDF(mergedDf, files[i])
        else:
            tmp = pd.read_table(files[i], low_memory=False, comment='#')
            tmp = AddColumnPercent(tmp, max(tmp['Cumulate_Sample_Number']), 'Covered_Sample_Number', 'Covered_Sample_Number_Percent')
            tmp.columns = RenameDF(tmp, files[i])
            # mode = union:outer or intersect:inner
            mergedDf = mergedDf.merge(tmp, on=cols, how=mode)    
    # first sort by covered_sample sum, then re-calculate cumulate_sample_sum and cumulate_sample_sum_percent
    covered_sample_cols = ColumnsByPattern(mergedDf.columns, 'Covered_Sample_Number')
    covered_sample_sum = mergedDf[covered_sample_cols].sum(axis=1, skipna=True, numeric_only=None)
    covered_sample_percent_cols = ColumnsByPattern(mergedDf.columns, 'Covered_Sample_Number_Percent')
    covered_sample_percent_sum = mergedDf[covered_sample_percent_cols].sum(axis=1, skipna=True, numeric_only=None)
    covered_sample_percent_average = covered_sample_percent_sum/len(covered_sample_percent_cols)
    mergedDf = pd.concat([mergedDf, pd.DataFrame({'Covered_Sample_Sum':covered_sample_sum, 'Average_Covered_Sample_Number_Percent':covered_sample_percent_average})], axis=1)
    if sort == 'CoveredSampleSum':
        mergedDf = mergedDf.sort_values('Covered_Sample_Number', ascending=False).reset_index().drop('index', axis=1).fillna('')
    elif sort == 'AverageCoveredSamplePercent':
        mergedDf = mergedDf.sort_values('Average_Covered_Sample_Number_Percent', ascending=False).reset_index().drop('index', axis=1).fillna('')
    # ------------------------- #
    # re-calculate cumulate sample and percent for each site
    merged_covered_sample = [[] for i in range(0,mergedDf.shape[0])]
    for f in files:
        covered_sample_col = f+':Covered_Sample'
        covered_sample_number_col = f+':Covered_Sample_Number'
        merged_sampleid = []
        merged_sampleid_numbers_cumulate = []
        for s in mergedDf.index:
            if mergedDf.iloc[s,][covered_sample_number_col] != '':
                # extract sites by index and compare with submaf
                # get sample ids
                covered_sample = mergedDf.iloc[s,][covered_sample_col]
                if type(covered_sample) is set:
                    merged_sampleid.append(covered_sample)
                    # merged sample id as new column 'Covered_Sample'
                    merged_covered_sample[s] = set(list(merged_covered_sample[s])+list(covered_sample))
                else:
                    # sample id is string type in 'set()' display format
                    s_sampleid = set([i.strip() for i in re.sub("['{}]","",covered_sample).split(',')])
                    merged_sampleid.append(s_sampleid)
                    # merged sample id as new column 'Covered_Sample'
                    merged_covered_sample[s] = set(list(merged_covered_sample[s])+list(s_sampleid))
                if merged_sampleid_numbers_cumulate == []:
                    merged_sampleid_numbers_cumulate.append(len(merged_sampleid[0]))
                else:
                    merged_sampleid_numbers_cumulate.append(len(set([j for i in merged_sampleid for j in i])))
            else:
                merged_sampleid.append('')
                if merged_sampleid_numbers_cumulate == []:
                    merged_sampleid_numbers_cumulate.append(0)
                else:
                    merged_sampleid_numbers_cumulate.append(merged_sampleid_numbers_cumulate[-1])
        mergedDf[f+":Cumulate_Sample_Number"] = merged_sampleid_numbers_cumulate
        mergedDf = AddColumnPercent(mergedDf, max(mergedDf[f+':Cumulate_Sample_Number']), f+':Cumulate_Sample_Number', f+':Cumulate_Sample_Number_Percent')
    # add merged sample ids
    mergedDf['Covered_Sample'] = merged_covered_sample
    mergedDf['Covered_Sample_Number'] = [len(i) for i in merged_covered_sample]
    # add multiple Cumulate_Sample_Number
    # Cumulate_Sample_Number
    cumulate_sample_cols = ColumnsByPattern(mergedDf.columns, 'Cumulate_Sample_Number')
    cumulate_sample_sum = mergedDf[cumulate_sample_cols].sum(axis=1, skipna=True, numeric_only=None)
    mergedDf = pd.concat([mergedDf, pd.DataFrame({'Cumulate_Sample_Number':cumulate_sample_sum})], axis=1)
    # Cumulate_Sample_Number_Percent
    mergedDf = AddColumnPercent(mergedDf, max(mergedDf['Cumulate_Sample_Number']), 'Cumulate_Sample_Number', 'Cumulate_Sample_Number_Percent')
    # output new dataframe
    mergedDf.to_csv(out, sep='\t', index=False, mode='w')


def main():
    args = GetArgs()
    if args.mode == 'union':
        mode = 'outer'
    elif args.mode == 'intersect':
        mode = 'inner'
    MergeDataframes(args.files, args.out, mode, args.sort)
    

if __name__ == '__main__':
    main()

