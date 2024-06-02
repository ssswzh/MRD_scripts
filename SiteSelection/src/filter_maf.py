#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2022/06/17
# @Author  : zhangsiwen
# @contact : zhang.siwen
# @LastModified: 2022/09/09
# @ChangeLog
#     20220617, first version
#     20220620, change optional argument to None, add ConditionFilter(), BedFilter()
#     20220629, separate othermax, correct position for indels in BedFilter()
#     20220909, add --func --indel, change --dbaf to --prev


import os
import argparse
import pandas as pd


def GetArgs():
    parser = argparse.ArgumentParser(description='Filter MAF-like file by arguments', formatter_class=argparse.RawTextHelpFormatter)
    parser._action_groups.pop()
    required = parser.add_argument_group('Required arguments')
    required.add_argument('--maf', help='maf-like file (dataframe with header), lines start with "#" will be igored', action='store', dest='maf', required=True)
    required.add_argument('--out', help='out file name', action='store', dest='out', required=True)
    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument('--bed', help='keep variants in bed file', default=None, action='store', dest='bed', required=False)
    # tumor
    optional.add_argument('--tdp', help='tumor total depth column name and cutoff (keep variant whose total depth is >= this), \n    separated by colon ":", default None (not filter), \n    e.g. "t_depth:100", "tumor_DP:100"', default=None, action='store', dest='tdp', required=False)
    optional.add_argument('--tad', help='tumor variant depth column name and cutoff (keep variant whose total depth is >= this), \n    separated by colon ":", default None (not filter), \n    e.g. "t_alt_count:3", "tumor_VD:3"', default=None, action='store', dest='tad', required=False)
    optional.add_argument('--taf', help='tumor variant frequency column name and cutoff (keep variant whose frequency is >= this), \n    if give a float number, will automately calculate AF firstly then filtering, \n    separated by colon ":", default None (not filter), \n    e.g. "0.01", "tumor_AF:0.01"', default=None, action='store', dest='taf', required=False)
    # normal
    optional.add_argument('--ndp', help='normal total depth column name and cutoff (keep variant whose variant depth is >= this), \n    separated by colon ":", default None (not filter), \n    e.g. "n_depth:100", "normal_DP:50"', default=None, action='store', dest='ndp', required=False)
    optional.add_argument('--nad', help='normal variant depth column name and cutoff (keep variant whose variant depth is <= this), \n    separated by colon ":", default None (not filter), \n    e.g. "n_alt_count:3", "normal_VD:2"', default=None, action='store', dest='nad', required=False)
    optional.add_argument('--naf', help='tumor variant frequency column name and cutoff (keep variant whose frequency is <= this), \n    if give a float number, will automately calculate AF firstly, \n    separated by colon ":", default None (not filter), \n    e.g. "0.01", "normal_AF:0.01"', default=None, action='store', dest='naf', required=False)
    # other
    optional.add_argument('--vt', help='variant type, multiple types separated by comma ",", choices ["All","Complex","Deletion","Insertion","SNV"] \n    default None (not filter)\n    e.g. for multiple types: "TYPE:SNV,Insertion,Deletion", "TYPE:All" ', default=None, action='store', dest='vt', required=False)
    optional.add_argument('--prev', help='variant prevalence column pattern and cutoff in gnomAD, ExAC, and other database (keep variant whose frequency is <= this), \n    default None (not filter)\n    e.g. "_AF:0.01"', default=None, action='store', dest='prev', required=False)
    optional.add_argument('--othermin', help='other column pattern and cutoff (keep variant whose value is >= this), \n    separated by semi-colon ";", default None (not filter)\n    e.g. "tumor_QUAL:20;normal_QUAL:20;tumor_MQ:10;normal_MQ:10"', default=None, action='store', dest='othermin', required=False)
    optional.add_argument('--othermax', help='other column pattern and cutoff (keep variant whose value is <= this), \n    separated by semi-colon ";", default None (not filter)\n    e.g. "tumor_NM:2;normal_NM:2"', default=None, action='store', dest='othermax', required=False)
    optional.add_argument('--indel', help='if give maf, filter indel with size larger than this, default None', default=None, type=int, action='store', dest='indel', required=False)
    optional.add_argument('--func', help='if give maf, the level to filter mutations functions, default [high] \nlow: Frame_Shift_Del,Frame_Shift_Ins,Missense_Mutation,Nonsense_Mutation,Nonstop_Mutation,Silent,Splice_Site,Translation_Start_Site, \nhigh: Missense_Mutation,Nonsense_Mutation,Nonstop_Mutation, \nUse "None" for not filtering.', default='high' ,action='store', dest='func', required=False)
    usage = '''Usage:
For standard maf file:
    %(prog)s --maf maf --out out --bed bed --tdp t_depth:100 --tad t_alt_count:3 --taf 0.01 --ndp n_depth:50 --nad n_alt_count:2 --naf 0.01 --vt TYPE:All --prev _AF:0.01 > log.stats
For non-standard file:
    %(prog)s --maf maf --out out --bed bed --tdp tumor_DP:100 --tad tumor_VD:3 --taf tumor_AF:0.01 --ndp normal_DP:50 --nad normal_VD:2 --naf normal_AF:0.01 --vt TYPE:All --prev _AF:0.01 --othermin "tumor_QUAL:20;normal_QUAL:20;tumor_MQ:10;normal_MQ:10" --othermax "tumor_NM:2;normal_NM:2" > log.stats
For filter BED only:
    %(prog)s --maf maf --out out --bed bed > log.stats
\n'''
    parser.epilog = usage
    args = parser.parse_args()
    return args


'''
maffile = 'P04.zhengu.xls'
outfile = 'test'
bedfile = '/mnt/ddngs/zhangsw/project/MRD/OV_WES/bed/xgen-exome-research-panel-targets.b37.bed'
tdp = 'tumor_DP:100'
tad = 'tumor_VD:3'
taf = 'tumor_AF:0.01'
ndp = 'normal_DP:50'
nad = 'normal_VD:2'
naf = 'normal_AF:0.01'
vt = 'TYPE:All'
prev = '_AF:0.01'
othermin = 'tumor_QUAL:20;normal_QUAL:20;tumor_MQ:10;normal_MQ:10'
othermax = 'tumor_NM:2;normal_NM:2'
'''


Variant_Classification_low = ["Frame_Shift_Del","Frame_Shift_Ins","Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation","Silent","Splice_Site","Translation_Start_Site"]
Variant_Classification_high = ["Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation"]
Consequence_low = ['coding_sequence_variant','frameshift_variant','missense_variant','protein_altering_variant','splice_acceptor_variant','splice_donor_variant','start_lost','stop_gained','stop_lost','stop_retained_variant','synonymous_variant','start_retained_variant']
Consequence_high = ['missense_variant','protein_altering_variant','start_lost','stop_gained','stop_lost']


def PrintFilterInfo(info, number):
    '''
    @ func: print filter condition and left record number
    '''
    print('Filter ' + info + '\t' + str(number))


def ConditionFilter(df, operator, conditions): # 20220620
    '''
    @ func: filter df according to given argument
    @ parameter:
        conditions: Dict{}
        operator: 'LargerThan' OR 'SmallerThan'
    '''
    choices = {'LargerThan':'>=', 'SmallerThan':'<='}
    for k,v in conditions.items():
        if operator == 'LargerThan':
            df = df[df[k] >= v]
        elif operator == 'SmallerThan':
            df = df[df[k] <= v]
        PrintFilterInfo(k+choices[operator]+str(v), df.shape[0])
    return df
    

def BedFilter(bedfile, df): # 20220620
    '''
    @ func: filter df according to bed region
    '''
    bed = pd.read_csv(bedfile, low_memory=False, sep='\t', comment='#', header=None)
    bed_cols = ['chr','start','end','gene','score','strand','startP','endP']
    bed.columns = bed_cols[0:bed.shape[1]]
    bed['chr'] = bed['chr'].str.replace(r'chr', '')
    # assign position column name
    pos = 'Start_Position' if 'Start_Position' in df.columns else 'Position'
    overlap_indices = []
    for idx in df.index:
        # for small deletions, correct position, 20220629
        if df.loc[idx,'Tumor_Seq_Allele2'] == '-':
            # 1:214637984:CCAGCCTCTGGGCCA>C  TO  1:214637985:CAGCCTCTGGGCCA>-
            subbed = bed[(bed['chr']==df.loc[idx,'Chromosome']) & (bed['start']+1<=df.loc[idx,pos]-1) & (bed['end']>=df.loc[idx,pos]-1)]
        else:
            subbed = bed[(bed['chr']==df.loc[idx,'Chromosome']) & (bed['start']+1<=df.loc[idx,pos]) & (bed['end']>=df.loc[idx,pos])]
        if not subbed.empty:
            overlap_indices.append(idx)
            df.loc[idx,'bed_chr'] = list(subbed['chr'])[0]
            df.loc[idx,'bed_start'] = list(subbed['start'])[0]
            df.loc[idx,'bed_end'] = list(subbed['end'])[0]
        else:
            df = df.drop(idx)
    df = df.astype({'bed_chr':'int', 'bed_start':'int', 'bed_end':'int'}, errors='ignore')
    return df


def MutationFunctionFilter(mafDf, functions): # 20220909
    '''
    @ func: if give MAF, filter mutation functions
    functions = [Variant_Classification, Consequence]
    '''
    # filter mutation classification
    if 'Variant_Classification' in mafDf.columns: #20220617
        mafDf = mafDf[mafDf['Variant_Classification'].isin(functions[0])]
    # non-standard maf file don't have Variant_Classification column, use Consequence column to filter for functions
    elif 'Consequence' in mafDf.columns: #20220617
        consequence_list = [c.split('&') for c in mafDf['Consequence']]
        consequence_result = []
        for cons in consequence_list:
            func = set(True for i in cons if i in functions[1])
            if True in func:
                consequence_result.append('selected')
            else:
                consequence_result.append('unselected')
        mafDf['Consequence_functions'] = consequence_result
        mafDf = mafDf[mafDf['Consequence_functions']=='selected']
        del mafDf['Consequence_functions']
    return mafDf


def FilterMaf(maffile, outfile, bedfile, tdp, tad, taf, ndp, nad, naf, vt, prev, othermin, othermax, indel, func):
    # -- read file --
    maf = pd.read_csv(maffile, low_memory=False, sep='\t', comment='#')
    maf['Chromosome'] = maf['Chromosome'].astype(str).str.replace(r'chr', '')
    PrintFilterInfo('None', maf.shape[0])
    # -- filter bed file --
    if bedfile != None:
        maf = BedFilter(bedfile, maf)
        PrintFilterInfo('BED '+os.path.basename(bedfile), maf.shape[0])
    # -- indel length -- 20220909
    if indel != None:
        for idx in maf.index:
            maf.loc[idx, 'indel_len'] = max(len(maf.loc[idx, 'Reference_Allele']), len(maf.loc[idx, 'Tumor_Seq_Allele2']))
        maf = maf[maf['indel_len'] <= indel]
        del maf['indel_len']
        PrintFilterInfo('indel_length<='+str(indel), maf.shape[0])
    # -- filter variant type --
    if vt != None:
        vt = {vt.split(':')[0]:vt.split(':')[1].split(',')}
        if 'All' not in list(vt.values())[0]:
            maf = maf[maf[list(vt.keys())[0]].isin(list(vt.values())[0])]
        PrintFilterInfo(str(list(vt.keys())[0])+' in '+'|'.join(list(vt.values())[0]), maf.shape[0])
    # -- filter dp, ad, af -- # 20220620
    for element in [tdp, tad, taf, ndp]:
        if element is not None:
            if element==taf and ':' not in element:
                maf['t_allele_frequency'] = maf[tad.split(':')[0]]/maf[tdp.split(':')[0]]
                argumentDict = {'t_allele_frequency':float(element)}
            else:
                argumentDict = {element.split(':')[0]:float(element.split(':')[1])}
            maf = ConditionFilter(maf, 'LargerThan', argumentDict)
    # -- other arguments, keep variant whose value is smaller than cutoff -- # 20220620
    for element in [nad, naf]:
        if element is not None:
            if element==naf and ':' not in element:
                maf['n_allele_frequency'] = maf[nad.split(':')[0]]/maf[ndp.split(':')[0]]
                adaf = {'n_allele_frequency':float(element)}
            else:
                adaf = {element.split(':')[0]:float(element.split(':')[1])}
            maf = ConditionFilter(maf, 'SmallerThan', adaf)
    # -- other arguments, keep variant whose value is smaller than cutoff -- # 20220629
    if othermax is not None:
        othermax = {d.split(':')[0]:float(d.split(':')[1]) for d in othermax.split(';')}
        maf = ConditionFilter(maf, 'SmallerThan', othermax)
    # -- other arguments, keep variant whose value is larger than cutoff -- # 20220620
    if othermin is not None:
        othermin = {d.split(':')[0]:float(d.split(':')[1]) for d in othermin.split(';')}
        maf = ConditionFilter(maf, 'LargerThan', othermin)
    # -- filter af in databases -- # 20220620
    if prev is not None:
        prev = prev.split(':')
        if taf is not None and naf is not None:
            prev_cols = list(set(col for col in maf.columns if prev[0] in col) - set(n.split(':')[0] for n in [taf, naf]))
        else:
            prev_cols = list(set(col for col in maf.columns if prev[0] in col))
        for c in prev_cols:
            maf = maf[(maf[c].astype(float) <= float(prev[1])) | (maf[c].isnull())]
        PrintFilterInfo('all database AF<='+str(prev[1]), maf.shape[0])
    # -- filter variant functions --
    if func == 'high':
        functions = [Variant_Classification_high, Consequence_high]
        maf = MutationFunctionFilter(maf, functions)
        PrintFilterInfo('Variant_Classification in high', maf.shape[0])
    elif func == 'low':
        functions = [Variant_Classification_low, Consequence_low]
        maf = MutationFunctionFilter(maf, functions)
        PrintFilterInfo('Variant_Classification in low', maf.shape[0])
    # -- write output --
    out = open(outfile,'w')
    out.write('#version 2.4\n')
    out.close()
    maf.to_csv(outfile, sep='\t', index=False, mode='a')
    

def main():
    ''' '''
    args = GetArgs()
    FilterMaf(args.maf, args.out, args.bed, args.tdp, args.tad, args.taf, args.ndp, args.nad, args.naf, args.vt, args.prev, args.othermin, args.othermax, args.indel, args.func)


if __name__ == '__main__':
    main()
