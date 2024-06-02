#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2022/06/02
# @Author  : zhangsiwen
# @contact : zhang.siwen
# @LastModified: 2022/08/04
# @ChangeLog
#     20220602, first version
#     20220615, change final output format
#     20220707, add depth-of-read uniformity (measured as the ratio of the 80th percentile to the 20th percentile), defined in TRACERx 2017 Nature Signatera doi:10.1038/nature22364
#     20220802, Target_Bases>=0.2xMean_Coverage_Pct to [0.1,1], step 0.1
#     20220804, add Min_Target_Reads & Max_Target_Reads
#     20221008, picard will auto merge overlapped regions, so change total region numbers to picard result rather than bed region provided


import os
import argparse
from collections import OrderedDict
import pandas as pd
import numpy as np


def GetArgs():
    parser = argparse.ArgumentParser(description='QC stats for AmpliSeq', formatter_class=argparse.RawTextHelpFormatter)
    parser._action_groups.pop()
    required = parser.add_argument_group('Required arguments')
    required.add_argument('--metrics', help='Picard CollectTargetedPcrMetrics result', action='store', dest='metrics', required=True)
    required.add_argument('--per_base_cov', help='PER_BASE_COVERAGE from Picard CollectTargetedPcrMetrics', action='store', dest='per_base_cov', required=True)
    required.add_argument('--per_target_cov', help='PER_TARGET_COVERAGE from Picard CollectTargetedPcrMetrics', action='store', dest='per_target_cov', required=True)
    required.add_argument('--bed', help='amplicon bed file', action='store', dest='bed', required=True)
    required.add_argument('--out', help='out file name', action='store', dest='out', required=True)
    #optional = parser.add_argument_group('Optional arguments')
    args = parser.parse_args()
    return args

'''
metrics_file = '120030268.CollectTargetedPcrMetrics'
target_cov_file = '120030268.per_target_cov'
base_cov_file = '120030268.per_base_cov'
bed_file = '/mnt/ddngs/zhangsw/project/MRD/ThermoS5/bed/Oncomine_Myeloid.20200429.designed.merged.bed'
outfile = 'test'
'''

# report reads and base cutoff
reads_cutoff = [1,2,10,50,100,500,1000,2500,5000,10000,25000,50000,100000]
base_coverage_cutoff = [1,2,10,50,100,500,1000,2500,5000,10000,25000,50000,100000]
uniformity_cutoff = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]


def TargetedPcrMetrics(metrics_file, bed_file):
    '''
    func: extract info from metrics file
    return: dict{} 
    '''
    amplicon_number = int(os.popen("grep -Ev '^track|^$' %s |wc -l" % bed_file).readline().strip("\n"))
    metrics = open(metrics_file)
    final_stats = OrderedDict()
    for line in metrics:
        if line.startswith('#') or line.strip()=='':
            continue
        elif line.startswith('CUSTOM_AMPLICON_SET'):
            attributes = line.strip().split('\t')
        elif line.startswith('coverage_or_base_quality'):
            break
        else:
            values = line.strip().split('\t')
            metrics_dict = dict(zip(attributes,values))
            # assay info
            final_stats['Custom_Assay'] = os.popen("basename -s .bed %s" % bed_file).readline().strip("\n")
            final_stats['Amplicon_Numbers'] = str(amplicon_number)
            final_stats['Amplicon_Total_Size'] = metrics_dict['TARGET_TERRITORY']
            final_stats['Mean_Amplicon_Length'] = str(int(metrics_dict['TARGET_TERRITORY'])/amplicon_number)
            #final_stats['Total_Reads'] = metrics_dict['TOTAL_READS']
            final_stats['PF_Reads'] = metrics_dict['PF_READS']
            final_stats['PF_Bases'] = metrics_dict['PF_BASES']
            final_stats['Mapped_Reads'] = metrics_dict['PF_UQ_READS_ALIGNED']
            final_stats['Mapped_Reads_Pct'] = str(float(int(final_stats['Mapped_Reads'])/int(final_stats['PF_Reads'])))
            final_stats['Mapped_Bases'] = metrics_dict['PF_UQ_BASES_ALIGNED']
            final_stats['Mapped_Bases_Pct'] = str(float(int(final_stats['Mapped_Bases'])/int(final_stats['PF_Bases'])))
            # target info
            final_stats['Target_Reads'] = ''
            final_stats['Target_Reads_Pct'] = ''
            final_stats['Target_Bases'] = metrics_dict['ON_AMPLICON_BASES']
            final_stats['Target_Bases_Pct'] = str(float(int(metrics_dict['ON_AMPLICON_BASES'])/int(final_stats['PF_Bases'])))
            final_stats['Mean_Target_Reads'] = ''
            # min median mean max coverage, 20220804
            final_stats['Min_Target_Coverage'] = metrics_dict['MIN_TARGET_COVERAGE']
            final_stats['Median_Target_Coverage'] = metrics_dict['MEDIAN_TARGET_COVERAGE']
            final_stats['Mean_Target_Coverage'] = metrics_dict['MEAN_TARGET_COVERAGE']
            final_stats['Max_Target_Coverage'] = metrics_dict['MAX_TARGET_COVERAGE']
            final_stats['Target_Reads(80/20percentile)'] = '' # 20220707
            #final_stats['Target_Bases>=0.2xMean_Coverage_Pct'] = ''
            final_stats['Non-zero_Coverage_Target_Pct'] = str(1-float(metrics_dict['ZERO_CVG_TARGETS_PCT']))
            # target reads and base
            for read_number in reads_cutoff:
                keyname = 'Target_Reads>=' + str(read_number) + '_Pct'
                final_stats[keyname] = ''
            for base_coverage in base_coverage_cutoff:
                final_keyname = 'Target_Bases>=' + str(base_coverage) + 'X_Pct'
                metrics_keyname = 'PCT_TARGET_BASES_' + str(base_coverage) + 'X'
                final_stats[final_keyname] = metrics_dict[metrics_keyname]
            for pct in uniformity_cutoff: # 20220802
                keyname = 'Target_Bases>=' + str(pct) + 'xMean_Coverage_Pct'
                final_stats[keyname] = ''
    metrics.close()
    return final_stats


def PerTargetCov(final_stats, target_cov_file, base_cov_file):
    target_cov_matrix = pd.read_csv(target_cov_file, sep='\t', header=0, low_memory=False)
    final_stats['Target_Reads'] = str(sum(target_cov_matrix['read_count']))
    final_stats['Target_Reads_Pct'] = str(float(int(final_stats['Target_Reads'])/int(final_stats['PF_Reads'])))
    final_stats['Mean_Target_Reads'] = str(int(final_stats['Target_Reads'])/int(target_cov_matrix.shape[0])) # 20221008
    # depth-of-read uniformity (measured as the ratio of the 80th percentile to the 20th percentile), 20220707
    final_stats['Target_Reads(80/20percentile)'] = str(np.nanpercentile(target_cov_matrix['read_count'],20)/np.nanpercentile(target_cov_matrix['read_count'],80))
    # count target with reads number larger than cutoff
    for read_number in reads_cutoff:
        keyname = 'Target_Reads>=' + str(read_number) + '_Pct'
        target_number = target_cov_matrix[target_cov_matrix['read_count'].fillna(0).astype('int')>=read_number].shape[0]
        final_stats[keyname] = str(float(int(target_number)/int(target_cov_matrix.shape[0]))) # 20221008
    base_cov_matrix = pd.read_csv(base_cov_file, sep='\t', header=0, low_memory=False)
    for pct in uniformity_cutoff: # 20220802
        keyname = 'Target_Bases>=' + str(pct) + 'xMean_Coverage_Pct'
        coverage_cutoff = pct * float(final_stats['Mean_Target_Coverage'])
        final_stats[keyname] = str(float(base_cov_matrix[base_cov_matrix['coverage'].fillna(0).astype('int')>=coverage_cutoff].shape[0] / int(final_stats['Amplicon_Total_Size'])))
    return final_stats


def FinalOutput(final_stats, outfile):
    # write result
    out = open(outfile, 'w')
    for k,v in final_stats.items(): # 20220615
        out.write(k + '\t' + v + '\n')
    #out.write('\t'.join(list(final_stats.keys())) + '\n')
    #out.write('\t'.join(list(final_stats.values())) + '\n')
    out.close()
    

def main():
    args = GetArgs()
    final_stats = TargetedPcrMetrics(args.metrics, args.bed)
    final_stats = PerTargetCov(final_stats, args.per_target_cov, args.per_base_cov)
    FinalOutput(final_stats, args.out)
    

if __name__ == '__main__':
    main()
