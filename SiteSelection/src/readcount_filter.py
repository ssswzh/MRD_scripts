#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2022/05/23
# @Author  : zhangsiwen
# @contact : zhang.siwen
# @LastModified: 2022/05/23
# @ChangeLog
#     20220523, first version


import argparse
from collections import OrderedDict
import pandas as pd


def GetArgs():
    parser = argparse.ArgumentParser(description='Filter readcount_extraction file by add comment', formatter_class=argparse.RawTextHelpFormatter)
    parser._action_groups.pop()
    required = parser.add_argument_group('Required arguments')
    required.add_argument('--in', help='readcount_extraction output with "--detail" set to be True', action='store', dest='input', required=True)
    required.add_argument('--out', help='out file name', action='store', dest='out', required=True)
    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument('--dp', help='filter depth, default [500]', default=500, type=int, action='store', dest='dp', required=False)
    optional.add_argument('--ad', help='filter allele depth, default [2]', default=2, type=int, action='store', dest='ad', required=False)
    args = parser.parse_args()
    return args

'''
in = '/mnt/ddngs/zhangsw/project/MRD/OV_WES/plasma2/readcount/P103.bam-readcount.plasma_signal'
outfile = 'test.tsv'
dp = 500
ad = 2
'''


def FilterAndComment(infile, outfile, DPcutoff=500, ADcutoff=2):
    '''
    infile:
    	chrs	pos	sample	Ref:Alt:Depth:A:C:G:T:N
    1:5935125	1	5935125	0.054988	G:T:1273:0:0:1203:70:0
    1:6635151	1	6635151	0.181944	T:C:3622:0:658:1:2963:0
    1:28261721	1	28261721	0.059701	C:A:1005:59:945:1:0:0
    1:42900954	1	42900954	0.072165	T:C:1261:0:91:0:1170:0
    '''
    out = open(outfile, 'w')
    with open(infile) as readcount_extraction:
        for line in readcount_extraction:
            if line.startswith('\tchrs'):
                out.write(line.strip() + '\tFilter\n')
                continue
            mutation, chrs, pos, sample, info = line.strip().split('\t')
            details = dict(zip(['ref','alt','depth','A','C','G','T','N'], info.split(':')))
            comment = ''
            if int(details['depth']) < DPcutoff:
                comment = comment + 'DPcutoff' + str(DPcutoff)
            if int(details[details['alt']]) < ADcutoff:
                comment = comment + 'ADcutoff' + str(ADcutoff)
            comment = 'PASS' if comment=='' else comment
            out.write(line.strip() + '\t' + comment + '\n')
    out.close()


def main():
    ''' '''
    args = GetArgs()
    FilterAndComment(args.input, args.out, args.dp, args.ad)


if __name__ == '__main__':
    main()
