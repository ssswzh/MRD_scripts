#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2022/02/15
# @Author  : zhangsiwen
# @contact : zhang.siwen@puruijizhun.com
# @LastModified: 2022/03/09
# @ChangeLog
#     20220215, first version
#     20220309, if CNV region not exists, at least set the same tumor content for sample
#     20220310, add --correct
#     20220321, re-organize functions, add --source
#     20220323, add --purity, re-organize CombineDfPyclone()


import argparse
import pandas as pd
from scipy.stats import expon


def GetArgs():
    parser = argparse.ArgumentParser(description='Generate PyClone input file by vcf and CNV and clinic data', formatter_class=argparse.RawTextHelpFormatter)
    parser._action_groups.pop()
    required = parser.add_argument_group('Required arguments')
    required.add_argument('--vcf', help='somatic mutation vcf file', action='store', dest='vcf', required=True)
    required.add_argument('--cnv', help='cnv file by ASCAT', action='store', dest='cnv', required=True)
    required.add_argument('--out', help='out file name', action='store', dest='out', required=True)
    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument('--source', help='software to generate cnv file, default [facets]', default='facets', choices=['facets','ascat'], action='store', dest='source', required=False)
    optional.add_argument('--gender', help='clinic file, default [female]', default='female', choices=['male','female'], action='store', dest='gender', required=False)
    optional.add_argument('--sampleid', help='sampleid, default [sample]', default='sample', action='store', dest='sampleid', required=False)
    optional.add_argument('--depth', help='filter variant with total depth larger than this, default [30]', default=30, type=int, action='store', dest='depth', required=False)
    optional.add_argument('--ad', help='filter variant with allele depth larger than this, default [5]', default=5, type=int, action='store', dest='ad', required=False)
    optional.add_argument('--af', help='filter variant with allele frequency larger than this, default [0.01]', default=0.01, type=float, action='store', dest='af', required=False)
    optional.add_argument('--purity', help='give tumor purity manually', default='', action='store', dest='purity', required=False)
    #optional.add_argument('--correct', help='choose to correct for observed mutation copy number, default no correct', action='store_true', dest='correct', required=False)
    # usage examples
    usage = '''Usage:
    %(prog)s --vcf vcf --cnv cnv --out outdir --gender male --depth 30 --ad 5 --af 0.05 --source facets
'''
    parser.epilog = usage
    args = parser.parse_args()
    return args


'''
vcffile = 'ascat/P04.somatic.filterBlacklist.norm.vcf'
cnvfile = 'ascat/P04.cnv.ascat'
gender = 'female'
out = 'test'
depth_cutoff = 30
ad_cutoff = 5
af_cutoff = 0.01
'''

# out file header: ['mutation_id','sample_id','ref_counts','alt_counts','normal_cn','major_cn','minor_cn','tumour_content']


def Vcf2Dataframe(vcffile, sampleid, depth_cutoff, ad_cutoff, af_cutoff):
    '''
    @ func: read vcf file and return a matrix
    @ result: DF(mutationDf).columns = ['mutation_id','sample_id','ref_counts','alt_counts','normal_cn','major_cn','minor_cn','tumour_content']
    '''
    mutationid = []
    chrom = []
    position = []
    ref_counts = []
    alt_counts = []
    vaf = []
    with open(vcffile) as vcf:
        for line in vcf:
            if line.startswith('#'):
                continue
            records = line.strip().split('\t')
            if len(records) == 11:
                CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, NORMAL, TUMOR = records[0:11]
                format_dict = dict(zip(FORMAT.split(':'), TUMOR.split(':')))
                # depth, allele depth, allele frequency not meet conditions
                dp = int(format_dict['DP'])
                if len(format_dict['AD'].split(',')) > 1:
                    rd = int(format_dict['AD'].split(',')[0])
                    ad = int(format_dict['AD'].split(',')[1])
                else:
                    ad = int(format_dict['AD'])
                    rd = int(format_dict['RD'])
            elif len(records) == 13:
                CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, NORMAL1, TUMOR1, NORMAL2, TUMOR2 = records[0:13]
                # varscan did not discover this mutation
                if TUMOR1.startswith('./.'):
                    format_dict = dict(zip(FORMAT.split(':'), TUMOR2.split(':')))
                    dp = int(format_dict['DP'])
                    rd = int(format_dict['AD'].split(',')[0])
                    ad = int(format_dict['AD'].split(',')[1])
                else:
                    format_dict = dict(zip(FORMAT.split(':'), TUMOR1.split(':')))
                    dp = int(format_dict['DP'])
                    ad = int(format_dict['AD'])
                    rd = int(format_dict['RD'])
            af = round(float(ad/dp),6)
            if not ( dp >= int(depth_cutoff) and ad >= int(ad_cutoff) and af >= float(af_cutoff) ):
                continue
            mutationid.append(':'.join([CHROM, POS, REF, ALT]))
            chrom.append(CHROM)
            position.append(POS)
            ref_counts.append(rd)
            alt_counts.append(ad)
            vaf.append(af)
    # mutationDf header: ['mutation_id','sample_id','ref_counts','alt_counts','normal_cn','major_cn','minor_cn','tumour_content']
    mutationDf = pd.DataFrame({'mutation_id':mutationid, 'CHROM':chrom, 'POS':position, 'sample_id':sampleid, 'ref_counts':ref_counts, 'alt_counts':alt_counts, 'vaf':vaf, 'major_cn':2, 'minor_cn':0, 'tumour_content':0.5}, index=mutationid)
    mutationDf = mutationDf.astype({'mutation_id':'str', 'CHROM':'str', 'POS':'int', 'sample_id':'str', 'ref_counts':'int', 'alt_counts':'int'}, errors='ignore')
    return mutationDf


def AssignCNbyGender(row, gender):
    '''
    @ func: return cn number by gender, invoked by CombineDfPyclone()
    '''
    if gender == 'male':
        if row['Chromosome'] in tuple(['X','chrX','Y','chrY','23','24']):
            return 1
        else:
            return 2
    else:
        return 2


def CNVheader(source):
    if source == 'ascat':
        #header = ['sample','chr','startpos','endpos','nMajor','nMinor','nAraw','nBraw','ploidy','ACF','goodnessOfFit','psi','cnTotal']
        return 'chr','startpos','endpos','nMajor','nMinor','ACF'
    elif source == 'facets':
        #header = ['chrom','seg','num.mark','nhet','cnlr.median','mafR','segclust','cnlr.median.clust','mafR.clust','start','end','cf.em	tcn.em','lcn.em']
        return 'chrom','start','end','mcn.em','lcn.em','purity'
    return


def CombineDfPyclone(mutationDf, cnvfile, out, gender, source='facets', purity=''):#, correct):
    '''
    @ func: merge dataframes (mutation and CNV) and output file format as PyClone suggested
    '''
    # cnv file
    cnvDf = pd.read_csv(cnvfile, low_memory=False, sep='\t', comment='#')
    chrom, start, end, major_cn, minor_cn, purityCol = CNVheader(source)
    cnvDf['distance'] = cnvDf[end] - cnvDf[start]
    # add normal_cn column, change sex chromosome to number
    mutationDf['normal_cn'] = mutationDf.apply(AssignCNbyGender, gender=gender, axis=1)
    # set the same tumor content for sample
    if purity == '':
        if not pd.isna(list(cnvDf[purityCol])[0]):
            purity = list(cnvDf[purityCol])[0]
        else:
            purity = 0.5
    mutationDf['tumour_content'] = purity
    # change sex chromosome to 23,24
    if 'X' in set(cnvDf[chrom]) or 'chrX' in set(cnvDf[chrom]):
        cnvDf[chrom] = cnvDf[chrom].replace(['X','Y'],['23','24']).replace(['chrX','chrY'],['23','24'])
    if 'X' in set(mutationDf['CHROM']) or 'chrX' in set(mutationDf['CHROM']):
        mutationDf['CHROM'] = mutationDf['CHROM'].replace(['X','Y'],['23','24']).replace(['chrX','chrY'],['23','24'])
    # type conversion
    cnvDf[chrom] = cnvDf[chrom].astype('int32')
    mutationDf['CHROM'] = mutationDf['CHROM'].astype('int32')
    # rename columns
    colnames = list(mutationDf.columns)
    if 't_ref_count' in colnames:
        colnames[colnames.index('t_ref_count')] = 'ref_counts'
    if 't_alt_count' in colnames:
        colnames[colnames.index('t_alt_count')] = 'alt_counts'
    mutationDf.columns = colnames
    # merge mutation and CNV dataframes
    for idx in mutationDf.index:
        subcnv = cnvDf[cnvDf[chrom]==mutationDf.loc[idx,'CHROM']]
        pos = mutationDf['POS'][idx]
        tmp = subcnv[(subcnv[start]<=pos) & (subcnv[end]>=pos)].sort_values(by='distance', ascending=True)
        if not tmp.empty:
            # in case major cn is zero
            if list(tmp[major_cn])[0] != 0 and not pd.isna(list(tmp[major_cn])[0]):
                mutationDf.loc[idx, 'major_cn'] = list(tmp[major_cn])[0]
            if not pd.isna(list(tmp[minor_cn])[0]):
                mutationDf.loc[idx, 'minor_cn'] = list(tmp[minor_cn])[0]
    mutationDf = mutationDf.sort_values(by=['mutation_id','CHROM','POS'], ascending=[True,True,True])
    del mutationDf['CHROM']
    del mutationDf['POS']
    '''
    # choose to correct for copy number
    if correct:
        # observed copy number Nmut = VAF/p*[p*CNt+CNn*(1-p)]
        #mutationDf['Nmut'] = mutationDf['vaf']/mutationDf['tumour_content'] * (mutationDf['tumour_content']*(mutationDf['major_cn']+mutationDf['minor_cn']) + (1-mutationDf['tumour_content'])*mutationDf['normal_cn'] )
        # expected copy number 
        #mutationDf['Nchr'] = [round(i) for i in mutationDf['Nmut']/mutationDf['vaf']]
        mutationDf['major_cn'] = 2
        mutationDf['minor_cn'] = 0
        mutationDf['tumour_content'] = 0.5
        #mutationDf['vaf'] = (mutationDf['major_cn']+mutationDf['minor_cn'])*mutationDf['tumour_content'] / (mutationDf['tumour_content']*mutationDf['major_cn'] + (1-mutationDf['tumour_content'])*mutationDf['minor_cn'] )    
        #mutationDf['ref_counts'] = round(2*mutationDf['alt_counts']/mutationDf['vaf']-mutationDf['alt_counts']
    '''
    mutationDf.to_csv(out, sep='\t', index=False, mode='w')


def main():
    ''' '''
    args = GetArgs()
    # mutationDf = Vcf2Dataframe(vcffile, 'sample', depth_cutoff, ad_cutoff, af_cutoff)
    mutationDf = Vcf2Dataframe(args.vcf, args.sampleid, args.depth, args.ad, args.af)
    CombineDfPyclone(mutationDf, args.cnv, args.out, args.gender, args.source, args.purity)#, args.correct)


if __name__ == '__main__':
    main()
