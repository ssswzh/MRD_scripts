#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2022/01/25
# @Author  : zhangsiwen
# @contact : zhang.siwen@puruijizhun.com
# @LastModified: 2022/01/25
# @ChangeLog
#     20220125, first version


import os
import argparse
import pysam
import datetime


def GetArgs():
    parser = argparse.ArgumentParser(description='Summarize tumor hotspot according to maf, compare hotspot according to baseline.', formatter_class=argparse.RawTextHelpFormatter)
    parser._action_groups.pop()
    required = parser.add_argument_group('Required arguments')
    required.add_argument('--vcf', help='vcf file', action='store', dest='vcf', required=True)
    required.add_argument('--out', help='output file', action='store', dest='out', required=True)
    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument('--depth', help='filter variant with total depth larger than this, default [30]', default=30, type=int, action='store', dest='depth', required=False)
    optional.add_argument('--ad', help='filter variant with allele depth larger than this, default [5]', default=5, type=int, action='store', dest='ad', required=False)
    optional.add_argument('--af', help='filter variant with allele frequency larger than this, default [0.05]', default=0.05, type=float, action='store', dest='af', required=False)
    # usage examples
    usage = '''Usage:
    %(prog)s --vcf file.vcf --out filtered.vcf
'''
    parser.epilog = usage
    args = parser.parse_args()
    return args

'''
vcffile = 'CRUK0034-R1.vcf'
outfile = 'test.vcf'
depth = 30
alleledepth = 5
allelfrequency = 0.05
'''

'''
vcf = pysam.VariantFile(vcffile)
header = vcf.header
header.formats.add("AF", number=1, type='String', description='Allele frequency.')
output = pysam.VariantFile(outfile, mode='w', header=header)

for record in vcf:
    sampleid = record.samples.keys()[0]
    depth = sum(i for i in record.samples[sampleid]['AD'])
    frequency = tuple(i/depth for i in record.samples[sampleid]['AD'][1:])
'''    


def FilterVcf(vcffile, outfile, depth, alleledepth, allelfrequency, command):
    vcf = open(vcffile)
    out = open(outfile, 'w')
    for record in vcf:
        if record.startswith('##'):
            out.write(record)
            continue
        if record.startswith('#CHROM'):
            out.write('##FORMAT=<ID=AF,Number=1,Type=Float,Description="Allele Frequency">\n')
            out.write('##Command: python ' + command + '; ' + str(datetime.datetime.today()) + '\n')
            out.write(record)
            continue
        ele = record.strip().split('\t')
        ele[8] = 'GT:AD:DP:AF'
        gt, ad, dp = ele[9].split(':')
        if ad == '.,.':
            continue
        ad = ad.split(',')
        dp = sum(int(i) for i in ad if i!='.')
        af = [str(float(i)/dp) for i in ad[1:]]
        if dp>=depth and int(ad[1])>=alleledepth and float(af[0])>=allelfrequency:
            ele[6] = 'PASS'
        else:
            ele[6] = 'Filtered'
        ad = ','.join(ad)
        af = ','.join(af)
        ele[9] = ':'.join([gt,ad,str(dp),af])
        out.write('\t'.join(ele)+'\n')
    out.close()


def main():
    args = GetArgs()
    # command line arguments
    command = [os.path.basename(__file__)]
    for k, v in dict(args._get_kwargs()).items():
        if v != None:
            command.append('--'+k)
            command.append(str(v))
    command = ' '.join(command)
    FilterVcf(args.vcf, args.out, args.depth, args.ad, args.af, command)
    
    
if __name__ == '__main__':
    main()
