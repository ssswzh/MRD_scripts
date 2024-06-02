#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2022/01/10
# @Author  : zhangsiwen
# @contact : zhang.siwen@puruijizhun.com
# @LastModified: 2022/01/10
# @ChangeLog
#     20220110, first version


import pandas as pd
import argparse


def GetArgs():
    parser = argparse.ArgumentParser(description='Merge infomation according to data sample barcode to maf.', formatter_class=argparse.RawTextHelpFormatter)
    parser._action_groups.pop()
    required = parser.add_argument_group('Required arguments')
    required.add_argument('--maf', help='maf file', action='store', dest='maf', required=True)
    required.add_argument('--info', help="TCGA baseline, containing several columns", action='store', dest='info', required=True)
    required.add_argument('--out', help='output file name', action='store', dest='out', required=True)
    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument('--col', help='merge dataframes by column name in info file', action='store', dest='col', default=['SAMPLE_BARCODE'], required=False)
    args = parser.parse_args()
    return args


'''
maffile = 'data/data_mutations.txt'
infofile = 'data_clinical_sample.info.txt'
out = 'test'
col = 'SAMPLE_BARCODE'
'''


def ReadDataframe(inputfile):
    df = pd.read_table(inputfile, low_memory=False, comment='#')
    return df


def MergeDataframe(maf, info, col, out):
    merged_df = pd.merge(maf, info, left_on="Tumor_Sample_Barcode", right_on=col)
    merged_df.to_csv(out, sep="\t", index=False)


def main():
    ''' '''
    args = GetArgs()
    # read dataframes
    info = ReadDataframe(args.info)
    # filter maf by parameters
    maf = ReadDataframe(args.maf)
    # merge dataframes by col
    MergeDataframe(maf, info, args.col, args.out)


if __name__ == '__main__':
    main()
