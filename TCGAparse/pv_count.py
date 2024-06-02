#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2022/01/06
# @Author  : zhangsiwen
# @contact : zhang.siwen@puruijizhun.com
# @LastModified: 2022/01/06
# @ChangeLog
#     20220106, first version
#     20220113, add maf processing


import argparse
import pysam
import itertools
import pandas as pd


def GetArgs():
    parser = argparse.ArgumentParser(description='Summarize Phased Variant (PV) in vcf file.', formatter_class=argparse.RawTextHelpFormatter)
    parser._action_groups.pop()
    required = parser.add_argument_group('Required arguments')
    required.add_argument('--input', help='vcf or maf to be process', action='store', dest='input', required=True)
    required.add_argument('--out', help='output prefix', action='store', dest='out', required=True)
    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument('--format', help='file format, vcf or maf, default [vcf]', action='store', dest='format', default='vcf', required=False)
    optional.add_argument('--size', help='mutations in which size to be considered as PV, default [170]', action='store', dest='size', default=170, type=int, required=False)
    args = parser.parse_args()
    return args


'''
vcffile = 'vcfs/DO52574.vcf.gz'
maffile = 'CRC/crc_apc_msk_impact_2020/data_mutation.add_clinical.txt'
out = 'testPV'
size = 170
'''


def CombinationNumber(total, n):
    '''
    @ func : NOT USED. get combinations for a list and return combination number
    '''
    l = []
    for i in itertools.combinations(range(total), n):
        l.append(i)
    return len(l)


def ReadVcf(vcffile):
    '''
    @ func : reaf vcf and return all records
    '''
    vcf = pysam.VariantFile(vcffile)
    recordCollection = []
    for record in vcf:
        recordCollection.append(record)
    return recordCollection


def VcfPVCounting(recordCollection, size):
    '''
    @ func : count PV for all sites in vcf, return [[PV combination], [PV combination]]
    '''
    pvCollection = []
    for i in range(0,len(recordCollection)):
        first = recordCollection[i]
        tmppv = [first]
        # this site locates in the same chrom with first site
        for j in range(i+1,len(recordCollection)):
            this = recordCollection[j]
            if first.chrom == this.chrom and this.start-first.start <= size:
                tmppv.append(this)
            else:
                if len(tmppv) >= 2:
                    pvCollection.append(tmppv)
                    break
    return pvCollection


def PVnumbers(pvlen):
    '''
    @ func : get pv2 and pv3 numbers by pvlenth list [2, 2, 2, 3, 3, 2, 3, 2, 2]
    '''
    pv2 = pvlen.count(2)
    numbers = set(pvlen) - {2}
    pv3 = 0
    for n in numbers:
        pv3 = pv3 + pvlen.count(n)
    pv2 = pv2 + pv3
    return pv2, pv3


def MafPVCounting(maffile, out, size):
    '''
    @ func : count PV for all sites and samples in maf, return [[PV combination], [PV combination]]
    '''
    maf = pd.read_table(maffile, low_memory=False, comment='#')
    maf = maf.sort_values(['Tumor_Sample_Barcode','Chromosome','Start_Position'], ascending = [True, True, True])
    sampleCollection = list(maf['Tumor_Sample_Barcode'].unique())
    pvCollection = {}
    for sample in sampleCollection:
        submaf = maf[maf['Tumor_Sample_Barcode']==sample].reset_index().drop('index', axis=1)
        tmpCollection = []
        for i in range(0, submaf.shape[0]):
            first = submaf.iloc[i,].to_dict()
            tmppv = [first]
            for j in range(i+1,submaf.shape[0]):
                this = submaf.iloc[j,].to_dict()
                if (first['Chromosome'] == this['Chromosome']) and (this['Start_Position']-first['Start_Position'] <= size) and (this['Start_Position'] != first['Start_Position']):
                    tmppv.append(this)
                else:
                    if len(tmppv) >= 2:
                        tmpCollection.append(tmppv)
                        break
        if tmpCollection != []:
            pvCollection[sample] = tmpCollection
    # structure: pvCollection['P-0000682-T01-IM3'] = [[record1,record2],[record3,record4,record5],[],[]]
    '''
    for i in pvCollection.keys():
        print(i, len(pvCollection[i]), [len(j) for j in pvCollection[i]])
    P-0000561-T01-IM3 1 [2]
    P-0000616-T01-IM3 1 [2]
    P-0000625-T01-IM3 1 [3]
    P-0000682-T01-IM3 3 [2, 2, 2]
    P-0000788-T01-IM3 4 [3, 2, 2, 2]
    P-0000862-T01-IM3 9 [2, 2, 2, 3, 3, 2, 3, 2, 2]
    '''
    return pvCollection


def SummarizePV(pvCollection, outfile, filetype):
    '''
    @ func : output PV stats and PV combinations
    name = outfile.split('/')[-1].split('.')[0]
    # PV combination site numbers
    pvlen = [len(i) for i in pvCollection]
    pv2, pv3 = PVnumbers(pvlen)
    # PV stats
    out = open(outfile+'.PVcount.tsv', 'w')
    out.write('File\tnumPVs (2x)\tnumPvs (3x)\n')
    out.write('%s\t%s\t%s\n' % (name, pv2, pv3))
    out.close()
    # PV combinations
    out = open(outfile+'.PVcombination.tsv', 'w')
    out.write('File\tPV1\tPV2\tPV3\tPV4\n')
    for comb in pvCollection:
        info = []
        for i in comb:
            info.append("_".join([i.chrom, str(i.start+1), i.ref, "|".join(list(i.alts))]))
        out.write('%s\t%s\n' % (name, '\t'.join(info)))
    out.close()
    @ func : output PV stats and PV combinations
    '''
    outstats = open(outfile+'.PVcount.tsv', 'w')
    outstats.write('Sample\tnumPVs (2x)\tnumPvs (3x)\n')
    outcomb = open(outfile+'.PVcombination.tsv', 'w')
    outcomb.write('File\tPV1\tPV2\tPV3\tPV4\n')
    if filetype == 'vcf':
        name = outfile.split('/')[-1].split('.')[0]
        # PV combination site numbers
        pvlen = [len(i) for i in pvCollection]
        pv2, pv3 = PVnumbers(pvlen)
        # PV stats
        outstats.write('%s\t%s\t%s\n' % (name, pv2, pv3))
        # PV combinations
        for comb in pvCollection:
            info = []
            for i in comb:
                info.append("_".join([i.chrom, str(i.start+1), i.ref, "|".join(list(i.alts))]))
            outcomb.write('%s\t%s\n' % (name, '\t'.join(info)))
    elif filetype == 'maf':
        # PV stats
        for i in pvCollection.keys():
            # PV combination site numbers
            pvlen = [len(j) for j in pvCollection[i]]
            pv2, pv3 = PVnumbers(pvlen)
            outstats.write('%s\t%s\t%s\n' % (i, pv2, pv3))
            # PV combinations
            for comb in pvCollection[i]:
                info = []
                for ele in comb:
                    info.append("_".join([ele['Chromosome'], str(ele['Start_Position']), ele['Reference_Allele'], '|'.join([ele['Tumor_Seq_Allele1'],ele['Tumor_Seq_Allele2']])]))
                outcomb.write('%s\t%s\t%s\t%s\n' % (i, ele['DISEASE'], ele['SUBTYPE'], '\t'.join(info)))
    outstats.close()
    outcomb.close()


def main():
    ''' '''
    args = GetArgs()
    if args.format == 'vcf':
        recordCollection = ReadVcf(args.input)
        pvCollection = VcfPVCounting(recordCollection, args.size)
    if args.format == 'maf':
        pvCollection = MafPVCounting(args.input, args.out, args.size)
    SummarizePV(pvCollection, args.out, args.format)
    

if __name__ == '__main__':
    main()
