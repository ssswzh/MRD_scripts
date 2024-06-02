#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2022/03/17
# @Author  : zhangsiwen
# @contact : zhang.siwen
# @LastModified: 2022/05/24
# @ChangeLog
#     20220318, first version, build based on somatic_baf_by_bamcount.py
#     20220324, minor change, fix bug when freq=0
#     20220325, minor change, fix bug that freq is wrote twice (same as previous one) and column becomes 6
#     20220524, if position not in readcount file, still give ref and alt info A:C:0:0:0:0:0:0 rather than 0:0:0:0:0:0:0:0


import argparse
from collections import OrderedDict
import pandas as pd


def GetArgs():
    parser = argparse.ArgumentParser(description='Calculate allele frequency according to readcount file', formatter_class=argparse.RawTextHelpFormatter)
    parser._action_groups.pop()
    required = parser.add_argument_group('Required arguments')
    required.add_argument('--vcf', help='vcf file for certain sites, only accept snps', action='store', dest='vcf', required=True)
    required.add_argument('--count', help='bam count result, either from bam-readcount or mpileup2readcounts', action='store', dest='count', required=True)
    required.add_argument('--out', help='out file name', action='store', dest='out', required=True)
    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument('--mode', help='calculate alt frequency [vaf] or b-allele frequency [baf], default [vaf]', default='vaf', choices=['vaf','baf'], action='store', dest='mode', required=False)
    optional.add_argument('--soft', help='software used to generate --count file [bam-readcount] or [mpileup2readcounts], default [bam-readcount]', default='bam-readcount', choices=['bam-readcount','mpileup2readcounts'], action='store', dest='soft', required=False)
    optional.add_argument('--detail', help='give detail count for Ref:Alt:Depth:A:C:G:T:N', action='store_true', dest='detail', required=False)
    args = parser.parse_args()
    return args

'''
vcffile = '/mnt/ddngs/zhangsw/project/MRD/OV_WES/tissue/P04/ascat/P04.somatic.filterBlacklist.norm.onlysnp.vcf'
countfile = '/mnt/ddngs/zhangsw/project/MRD/OV_WES/tissue/P04/P04.somatic.bam-readcount'
outfile = 'test.tsv'
mode = 'vaf'
soft = 'bam-readcount'
detail = True
'''

bam_count_info = ['base',
'count',
'avg_mapping_quality',
'avg_basequality',
'avg_se_mapping_quality',
'num_plus_strand',
'num_minus_strand',
'avg_pos_as_fraction',
'avg_num_mismatches_as_fraction',
'avg_sum_mismatch_qualities',
'num_q2_containing_reads',
'avg_distance_to_q2_start_in_q2_reads',
'avg_clipped_length',
'avg_distance_to_effective_3p_end']

bam_count_field = ['chr',
'position',
'reference_base',
'depth',
bam_count_info
]

'''
mpileup2readcounts_field = ['chr',
'pos',
'depth',
'ref_base',
'refcount',
'altcount',
'acount',
'ccount',
'gcount',
'tcount',
'ncount',
'indelcount',
'identifier'
]
'''

def BamReadcount(countfile, variant_sites, variant_base, mode='vaf', detail=False):
    '''
    bam-readcount result:
        chr position  reference_base  depth  bam_count_info
        1   12080     G               0      =:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00  A:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00  C:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00  G:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00  T:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00  N:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00
    '''
    base_order = ['A','C','G','T','N']
    with open(countfile) as count:
        for line in count:
            record = line.strip().split('\t')
            chrom, position, reference_base, depth = record[0:4]
            depth = int(depth)
            # store base count for ['A','C','G','T','N']
            base_detail = ['0'] * len(base_order)
            if chrom in variant_sites.keys() and position in variant_sites[chrom].keys():
                # alt base in vcf file
                alt = variant_base[chrom][position]
                if depth != 0:
                    alt_base_sum = 0
                    #print(chrom, position, variant_sites[chrom][position], alt)
                    for base_data_string in record[5:9]:
                        base_values = base_data_string.split(':')
                        base_info = dict(zip(bam_count_info, base_values))
                        base_detail[base_order.index(base_info['base'])] = base_info['count']
                        if mode == 'baf':
                            if base_info['base'] != reference_base:
                                alt_base_sum = alt_base_sum + int(base_info['count'])
                            freq = str(round(float(alt_base_sum/depth),6))
                        elif mode == 'vaf':
                            if base_info['base'] == alt:
                                freq = str(round(float(int(base_info['count'])/depth),6))
                else:
                    depth = '0'
                    freq = '0'
                    base_detail = '0:0:0:0:0'
                if detail:
                    freq = freq + '\t' + ':'.join([reference_base,alt,str(depth)]) + ':' + ':'.join(base_detail)
                variant_sites[chrom][position] = freq
            else:
                continue
    return variant_sites


def MpileupReadcount(countfile, variant_sites, variant_base, mode='vaf', detail=False):
    '''
    mpileup2readcounts result:
        chr pos       depth  ref_base  refcount  altcount  acount  ccount  gcount  tcount  ncount  indelcount  identifier
        1   62349898  670    T         670       0         0       0       0       670     0       0           P32
    '''
    with open(countfile) as count:
        for line in count:
            record = line.strip().split('\t')
            chrom, position, depth, ref = record[0:4]
            if chrom == 'chr' and position == 'pos':
                mpileup2readcounts_field = record
                continue
            depth = int(depth)
            if chrom in variant_sites.keys() and position in variant_sites[chrom].keys():
                alt = variant_base[chrom][position]
                if depth != 0:
                    base_info = dict(zip(mpileup2readcounts_field[4:11], record[4:11]))
                    base_detail = ':'.join([base_info[i] for i in ['acount','ccount','gcount','tcount','ncount']])
                    if mode == 'baf':
                        freq = str(round(float(int(base_info['altcount'])/depth),6))
                        '''
                        alt_base_sum = 0
                        for k,v in base_info.items():
                            if k.upper() != reference_base:
                                alt_base_sum = alt_base_sum + int(base_info[k])
                        variant_sites[chrom][position] = str(round(float(alt_base_sum/depth),6))
                        '''
                    elif mode == 'vaf':
                        freq = str(round(float(int(base_info[alt.lower()+'count'])/depth),6))
                else:
                    depth = '0'
                    freq = '0'
                    base_detail = '0:0:0:0:0'
                if detail:
                    freq = freq + '\t' + ':'.join([ref,alt,str(depth)]) + ':' + base_detail
                variant_sites[chrom][position] = freq
            else:
                continue
    return variant_sites


def ReadCountExtraction(countfile, outfile, variant_sites, variant_base, soft='bam-readcount', mode='vaf', detail=False):
    # read bam-count and calculate frequencies
    if soft == 'bam-readcount':
        variant_sites = BamReadcount(countfile, variant_sites, variant_base, mode, detail)
    elif soft == 'mpileup2readcounts':
        variant_sites = MpileupReadcount(countfile, variant_sites, variant_base, mode, detail)
    # write result
    with open(outfile, 'w') as out:
        if detail:
            out.write('\tchrs\tpos\tsample\tRef:Alt:Depth:A:C:G:T:N\n')
        else:
            out.write('\tchrs\tpos\tsample\n')
        for chrom,info in variant_sites.items():
            for pos,freq in info.items():
                mutation_id = chrom + ':' + pos
                out.write('\t'.join([mutation_id, chrom, pos, freq]) + '\n')
    return


def VariantSite(vcffile, detail):
    '''
    @func: read, store variant sites, return dicts
    @result:
        variant_sites: OrderedDict([('1', {'1': '0', '2': '0'}, '2':{}, ...)])
        variant_sites: OrderedDict([('1', {'1': 'A', '2': 'C'}, '2':{}, ...)])
    '''
    # read and store variant sites
    variant_sites = OrderedDict()
    variant_base = OrderedDict()
    with open(vcffile) as vcf:
        for line in vcf:
            if line.startswith('#'):
                continue
            record = line.strip().split('\t')
            chrs, pos = record[0:2]
            if chrs not in variant_sites.keys():
                variant_sites[chrs] = {}
                variant_base[chrs] = {}
            if pos not in variant_sites[chrs]:
                if detail: # added 20220524
                    variant_sites[chrs][pos] = '0.0\t' + record[3] + ':' + record[4] + ':0:0:0:0:0:0'
                else:
                    variant_sites[chrs][pos] = '0.0'
                variant_base[chrs][pos] = record[4]
    return variant_sites, variant_base


def main():
    ''' '''
    args = GetArgs()
    variant_sites, variant_base = VariantSite(args.vcf, args.detail)
    bamcount_sites = ReadCountExtraction(args.count, args.out, variant_sites, variant_base, args.soft, args.mode, args.detail)


if __name__ == '__main__':
    main()
