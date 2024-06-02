#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2022/03/17
# @Author  : zhangsiwen
# @contact : zhang.siwen
# @LastModified: 2022/09/27
# @ChangeLog
#     20220318, first version, build based on somatic_baf_by_bamcount.py
#     20220324, minor change, fix bug when freq=0
#     20220325, minor change, fix bug that freq is wrote twice (same as previous one) and column becomes 6
#     20220524, if position not in readcount file, still give ref and alt info A:C:0:0:0:0:0:0 rather than 0:0:0:0:0:0:0:0
#     20220927, rebuild script, also add --strand to separate forward and reverse strand read count, change output columns, deprecated --detail
#     20221009, change all ref base to uppercase


import argparse
from collections import OrderedDict
import pandas as pd
from pyfaidx import Fasta


def GetArgs():
    parser = argparse.ArgumentParser(description='Calculate allele frequency according to readcount file', formatter_class=argparse.RawTextHelpFormatter)
    parser._action_groups.pop()
    required = parser.add_argument_group('Required arguments')
    required.add_argument('--count', help='bam count result, either from bam-readcount or mpileup2readcounts', action='store', dest='count', required=True)
    required.add_argument('--out', help='out file name', action='store', dest='out', required=True)
    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument('--range', help='bed or vcfs are accepted.\nif give vcf file, only accept snps', action='store', dest='range', required=False)
    optional.add_argument('--ref', help='if give bed file, should give reference fasta file, default %(default)s', default='/mnt/user/errand/pipeline/MC3/ref/Homo_sapiens_assembly19.fasta', action='store', dest='ref', required=False)
    optional.add_argument('--mode', help='for vcf file, calculate alt frequency [vaf] or b-allele frequency [baf], default %(default)s', default='vaf', choices=['vaf','baf'], action='store', dest='mode', required=False)
    optional.add_argument('--soft', help='software used to generate --count file [bam-readcount] or [mpileup2readcounts], default %(default)s', default='bam-readcount', choices=['bam-readcount','mpileup2readcounts'], action='store', dest='soft', required=False)
    #optional.add_argument('--detail', help='give detail count for Ref:Alt:Depth:A:C:G:T:N', action='store_true', dest='detail', required=False)
    optional.add_argument('--strand', help='for forward and reverse strand, give two columns of read counts, only useful for --soft=bam-readcount', action='store_true', dest='strand', required=False)
    usage = '''Usage:
  Get all positions:
    %(prog)s --count bam.readcount --out outfile
  Get positions inside bed file:
    %(prog)s --count bam.readcount --out outfile --range bed [--strand]
  Get positions inside vcf file:
    %(prog)s --count bam.readcount --out outfile --range vcf [--mode vaf/baf] [--strand]
    
'''
    parser.epilog = usage
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


def checkdepth(depth, ad):
    if depth == 0:
        return '0'
    else:
        return str(round(float(ad/depth),6))


def BamReadcount(countfile, sites, mode='vaf', strand=False):
    '''
    bam-readcount result:
        chr position  reference_base  depth  bam_count_info
        1	139409	A	10300	=:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00	A:10299:39.19:36.43:0.02:5496:4803:0.47:0.00:5.38:5496:0.52:141.57:0.45	C:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00	G:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00	T:1:28.00:37.00:0.00:0:1:0.99:0.06:114.00:0:0.00:142.00:0.47	N:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00
    '''
    base_order = ['A','C','G','T','N']
    result = sites.copy()
    count = open(countfile)
    for line in count:
        record = line.strip().split('\t')
        chrom, position, ref, depth = record[0:4]
        ref = ref.upper()
        depth = int(depth)
        # store base count for ['A','C','G','T','N'] # both forward and reverse 20220927
        base_detail_forward = [0] * len(base_order)
        base_detail_reverse = [0] * len(base_order)
        if depth != 0:
            for base_data_string in record[5:9]:
                base_values = base_data_string.split(':')
                base_info = dict(zip(bam_count_info, base_values))
                base_detail_forward[base_order.index(base_info['base'])] = int(base_info['num_plus_strand'])
                base_detail_reverse[base_order.index(base_info['base'])] = int(base_info['num_minus_strand'])
            depth_forward = sum(base_detail_forward)
            depth_reverse = sum(base_detail_reverse)   
            alt_base_sum_forward = depth_forward - int(base_detail_forward[base_order.index(ref)])
            alt_base_sum_reverse = depth_reverse - int(base_detail_reverse[base_order.index(ref)])
            total_base_detail = [i+j for i,j in zip(base_detail_forward,base_detail_reverse)]
            if sites: # region provided
                if chrom in sites.keys() and position in sites[chrom].keys():
                    # alt base in vcf file
                    alt = sites[chrom][position][1]
                    if alt != '?' and mode == 'vaf':
                        alt_f = int(base_detail_forward[base_order.index(alt)])
                        alt_r = int(base_detail_reverse[base_order.index(alt)])
                        freq_forward = checkdepth(depth_forward, alt_f)
                        freq_reverse = checkdepth(depth_reverse, alt_r)
                        freq = checkdepth(depth_forward, alt_f+alt_r)
                    else:
                        freq_forward = checkdepth(depth_forward, alt_base_sum_forward)
                        freq_reverse = checkdepth(depth_reverse, alt_base_sum_reverse)
                        freq = checkdepth(depth, alt_base_sum_forward+alt_base_sum_reverse)
                else:
                    continue
            else: # region NOT provided
                alt = '?'
                freq_forward = checkdepth(depth_forward, alt_base_sum_forward)
                freq_reverse = checkdepth(depth_reverse, alt_base_sum_reverse)
                freq = checkdepth(depth, alt_base_sum_forward+alt_base_sum_reverse)
            total = ':'.join([ref,alt,str(depth)]) + ':' + ':'.join([str(i) for i in total_base_detail])
            if strand:
                forward = ':'.join([ref,alt,str(depth_forward)]) + ':' + ':'.join([str(i) for i in base_detail_forward])
                reverse = ':'.join([ref,alt,str(depth_reverse)]) + ':' + ':'.join([str(i) for i in base_detail_reverse])
                freq = freq + '\t' + total + '\t' + '\t'.join([freq_forward, forward, freq_reverse, reverse])
            else:
                freq = freq + '\t' + total
        else:
            depth = '0'
            freq = '0'
            if strand:
                base_detail = ['0'] * 5
                info = ':'.join([ref,alt,str(depth)]) + ':' + ':'.join(base_detail)
                freq = freq + '\t' + info + '\t' + '\t'.join([freq, info, freq, info])
            else:
                freq = freq + '\t' + info
            result[chrom][position] = freq
        if chrom not in result.keys():
            result[chrom] = {}
        if position not in result[chrom]:
            result[chrom][position] = {}
        result[chrom][position] = freq
    count.close()
    return result, len(freq.split('\t'))


def MpileupReadcount(countfile, sites, mode='vaf'):
    '''
    mpileup2readcounts result:
        chr pos       depth  ref_base  refcount  altcount  acount  ccount  gcount  tcount  ncount  indelcount  identifier
        1   62349898  670    T         670       0         0       0       0       670     0       0           P32
    '''
    result = sites.copy()
    count = open(countfile)
    for line in count:
        record = line.strip().split('\t')
        chrom, position, depth, ref = record[0:4]
        if chrom not in result.keys():
            result[chrom] = {}
        if position not in result[chrom]:
            result[chrom][position] = {}
        if chrom == 'chr' and position == 'pos':
            mpileup2readcounts_field = record
            continue
        depth = int(depth)
        base_info = dict(zip(mpileup2readcounts_field[4:11], record[4:11]))
        base_detail = ':'.join([base_info[i] for i in ['acount','ccount','gcount','tcount','ncount']])
        if sites: #20220927
            if chrom in sites.keys() and position in sites[chrom].keys():
                alt = sites[chrom][position][1]
                if depth != 0:
                    if alt != '?' and mode == 'vaf':
                        freq = checkdepth(depth, base_info[alt.lower()+'count'])
                    else:
                        freq = checkdepth(depth, base_info['altcount'])
                else:
                    depth = '0'
                    freq = '0'
                    base_detail = '0:0:0:0:0'
                freq = freq + '\t' + ':'.join([ref,alt,str(depth)]) + ':' + base_detail
            else:
                continue
        else:
            alt = '?'
            freq = checkdepth(depth, base_info['altcount'])
            freq = freq + '\t' + ':'.join([ref,alt,str(depth)]) + ':' + base_detail
        result[chrom][position] = freq
    count.close()
    return result, len(freq.split('\t'))


def VariantSite(vcffile):
    '''
    @func: read, store variant sites, return dicts
    @result:
        sites: OrderedDict([('1', {'1': ['A','C'], '2': ['T','A']}, '2':{}, ...)])
    '''
    # read and store variant sites
    sites = OrderedDict()
    with open(vcffile) as vcf:
        for line in vcf:
            if line.startswith('#'):
                continue
            record = line.strip().split('\t')
            chrs, pos = record[0:2]
            if chrs not in sites.keys():
                sites[chrs] = {}
            if pos not in sites[chrs]: # 20220927
                sites[chrs][pos] = [record[3], record[4]]
    return sites


def BedSite(bedfile, reference):
    ''' 20220927
    @func: read, store variant sites, return dicts
    @result:
        sites: OrderedDict([('1', {'1': ['A','?'], '2': ['T','?']}, '2':{}, ...)])
    '''
    fasta = Fasta(reference, rebuild=False)
    sites = OrderedDict()
    with open(bedfile) as bed:
        for line in bed:
            if line.startswith('#'):
                continue
            record = line.strip().split('\t')
            chrs, start, end = record[0:3]
            if chrs not in sites.keys():
                sites[chrs] = {}
            for pos in range(int(start)+1,int(end)+1):
                ref = fasta[chrs][pos : pos+1].seq
                sites[chrs][str(pos)] = [ref, '?']
    return sites


def main():
    ''' '''
    args = GetArgs()
    # dependency
    if args.soft == 'mpileup2readcounts' and args.strand:
        print('Warning: mpileup2readcounts do not provide strand read count, ignore --strand and output only total read count.')
    if args.range: # 20220927
        if args.range.endswith('vcf'):
            sites = VariantSite(args.range)
        elif args.range.endswith('bed'):
            print('Warning: reference file is ' + args.ref)
            sites = BedSite(args.range, args.ref)
        else:
            print('Warning: not recognizing --range file type, ignoring file')
            sites = OrderedDict()
    else:
        sites = OrderedDict()
    
    # read bam-count and calculate frequencies
    if args.soft == 'bam-readcount':
        sites, col_num = BamReadcount(args.count, sites, args.mode, args.strand)
    elif args.soft == 'mpileup2readcounts':
        sites, col_num = MpileupReadcount(args.count, sites, args.mode)
    
    # write result
    with open(args.out, 'w') as out:
        if col_num == 2:
            header = ['','chrs','pos','af','Ref:Alt:Depth:A:C:G:T:N']
            out.write('\t'.join(header) + '\n')
        else:
            header = ['','chrs','pos','af','TotalRef:Alt:Depth:A:C:G:T:N', 'plus_af', 'PlusRef:Alt:Depth:A:C:G:T:N', 'minus_af', 'MinusRef:Alt:Depth:A:C:G:T:N']
            out.write('\t'.join(header) + '\n')
        for chrom,info in sites.items():
            for pos,freq in info.items():
                mutation_id = chrom + ':' + pos
                out.write('\t'.join([mutation_id, chrom, pos, freq]) + '\n')
    

if __name__ == '__main__':
    main()
