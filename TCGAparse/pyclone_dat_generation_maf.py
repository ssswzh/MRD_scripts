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
import pandas as pd


def GetArgs():
    parser = argparse.ArgumentParser(description='Generate PyClone input file by maf and CNV and clinic data', formatter_class=argparse.RawTextHelpFormatter)
    parser._action_groups.pop()
    required = parser.add_argument_group('Required arguments')
    required.add_argument('--maf', help='maf file', action='store', dest='maf', required=True)
    required.add_argument('--cnv', help='cnv file by ASCAT', action='store', dest='cnv', required=True)
    required.add_argument('--clin', help='clinic file', action='store', dest='clin', required=True)
    required.add_argument('--out', help='outdir', action='store', dest='out', required=True)
    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument('--by', help='output tsv by [Tumor_Sample_Barcode,PATIENT_BARCODE], default PATIENT_BARCODE', default='PATIENT_BARCODE', choices=['Tumor_Sample_Barcode','PATIENT_BARCODE'], action='store', dest='by', required=False)
    optional.add_argument('--depth', help='filter variant with total depth larger than this, default [30]', default=30, type=int, action='store', dest='depth', required=False)
    optional.add_argument('--ad', help='filter variant with allele depth larger than this, default [5]', default=5, type=int, action='store', dest='ad', required=False)
    optional.add_argument('--af', help='filter variant with allele frequency larger than this, default [0.05]', default=0.05, type=float, action='store', dest='af', required=False)
    optional.add_argument('--func', help="Variant_Classification type in maf file to be reported, separate by comma ',', use 'ALL' for NOT FILTER any classification \ndefault 'Frame_Shift_Del,Frame_Shift_Ins,In_Frame_Del,In_Frame_Ins,Intron,Missense_Mutation,Nonsense_Mutation,Nonstop_Mutation,RNA,Silent,Splice_Site,Translation_Start_Site' \nchoose from [ALL,3'Flank,3'UTR,5'Flank,5'UTR,Frame_Shift_Del,Frame_Shift_Ins,In_Frame_Del,In_Frame_Ins,Intron,Missense_Mutation,Nonsense_Mutation,Nonstop_Mutation,RNA,Silent,Splice_Site,Translation_Start_Site]", choices=['ALL','3\'Flank','3\'UTR','5\'Flank','5\'UTR','Frame_Shift_Del','Frame_Shift_Ins','In_Frame_Del','In_Frame_Ins','Intron','Missense_Mutation','Nonsense_Mutation','Nonstop_Mutation','RNA','Silent','Splice_Site','Translation_Start_Site'], action='store', dest='func', default="Frame_Shift_Del,Frame_Shift_Ins,In_Frame_Del,In_Frame_Ins,Intron,Missense_Mutation,Nonsense_Mutation,Nonstop_Mutation,RNA,Silent,Splice_Site,Translation_Start_Site", required=False)
    # usage examples
    usage = '''Usage:
    %(prog)s --maf maf --out outdir --by PATIENT_BARCODE --depth 30 --ad 5 --af 0.05
'''
    parser.epilog = usage
    args = parser.parse_args()
    return args


'''
maffile = 'data/data_mutation.add_clinical.txt'
cnvfile = 'data/data_copy_number.txt'
clinfile = 'data/data_clinical_patient.txt'
outdir = 'test'
by = 'Tumor_Sample_Barcode'
depth = 30
ad = 5
af = 0.05
'''

cols = ['mutation_id','sample_id','ref_counts','var_counts','vaf','major_cn','minor_cn','normal_cn','tumour_content','Variant_Classification','ObsPyCloneCCF','ObsPyCloneClonality','ObsPyCloneCluster','PhyloCCF','PyClonePhyloCCF','PyClonePhyloClonal','PyClonePhyloCluster','RegionSum']


def AssignCNbyGender(row, genderDict):
    '''
    @ func: return cn number by gender, invoked by CombineDfPyclone()
    '''
    if genderDict[row['PATIENT_BARCODE']] == 'Male':
        if row['Chromosome'] in tuple(['X','chrX','Y','chrY']):
            return 1
        else:
            return 2
    else:
        return 2


def CombineDfPyclone(maffile, clinfile, cnvfile, outdir, by, depth, ad, af):
    '''
    @ func: merge dataframes and output file format as PyClone suggested
    '''
    # if output dir exists
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    # read dataframes
    maf = pd.read_csv(maffile, low_memory=False, sep='\t', comment='#')
    clin = pd.read_csv(clinfile, low_memory=False, sep='\t', comment='#')
    gender = dict(zip(clin['PATIENT_ID'], clin['SEX']))
    # cnv file
    cnv = pd.read_csv(cnvfile, low_memory=False, sep='\t', comment='#')
    cnv['distance'] = cnv['endpos'] - cnv['startpos']
    if 'X' in set(cnv['chr']):
        cnv['chr'] = cnv['chr'].replace(['X','Y'],[23,24])
    # output file names
    outid = maf[by].unique()
    # add new columns
    maf = pd.concat([maf, pd.DataFrame({'major_cn':2,'minor_cn':0,'tumour_content':0.5}, index=maf.index)], axis=1)
    maf['normal_cn'] = maf.apply(AssignCNbyGender, genderDict=gender, axis=1)
    # add mutation_id and change column values (23 to X, 24 to Y)
    maf['mutation_id'] = maf[['Chromosome','Start_Position','Reference_Allele','Tumor_Seq_Allele2','Hugo_Symbol','HGVSc']].astype(str).agg(':'.join, axis=1)
    if 'X' in set(maf['Chromosome']):
        maf['Chromosome'] = maf['Chromosome'].replace(['X','Y'],[23,24])
    # rename columns
    colnames = list(maf.columns)
    if 't_ref_count' in colnames:
        colnames[colnames.index('t_ref_count')] = 'ref_counts'
    if 't_alt_count' in colnames:
        colnames[colnames.index('t_alt_count')] = 'var_counts'
    maf.columns = colnames
    # for each sample 
    for i in outid:
        submaf = maf[maf[by]==i]
        if submaf.empty:
            continue
        sampleid = submaf['SAMPLE_BARCODE'].unique()
        subcnv = cnv[cnv['sample'].isin(sampleid)]
        patientid = str(submaf['PATIENT_BARCODE'].unique()[0])
        # filter by depth, ad, af
        submaf['Depth'] = list(submaf['ref_counts'] + submaf['var_counts'])
        submaf['vaf'] = list(submaf['var_counts']/submaf['Depth'])
        if not patientid.endswith('ctDNA'):
            if '-' in submaf['Reference_Allele'] or '-' in submaf['Tumor_Seq_Allele2']:
                submaf = submaf[submaf['Depth']>50]
                submaf = submaf[submaf['var_counts']>10]
                submaf = submaf[submaf['vaf']>af]
            else:
                submaf = submaf[submaf['Depth']>depth]
                submaf = submaf[submaf['var_counts']>ad]
                submaf = submaf[submaf['vaf']>af]
        # for each site
        for idx in submaf.index:
            subsubcnv = subcnv[subcnv['chr']==int(submaf.loc[idx,'Chromosome'])]
            pos = submaf['Start_Position'][idx]
            tmp = subsubcnv[(subsubcnv['startpos']<=pos) & (subsubcnv['endpos']>=pos)].sort_values(by='distance', ascending=True)
            if not tmp.empty:
                # in case major cn is zero
                if list(tmp['nMajor'])[0] != 0:
                    submaf.loc[idx, 'major_cn'] = list(tmp['nMajor'])[0]
                submaf.loc[idx, 'minor_cn'] = list(tmp['nMinor'])[0]
                submaf.loc[idx, 'tumour_content'] = list(tmp['ACF'])[0]
            elif i in set(subcnv['sample']):
                submaf.loc[idx, 'tumour_content'] = list(subcnv[subcnv['sample']==i]['ACF'].unique())[0]
        submaf['sample_id'] = list(submaf['Tumor_Sample_Barcode'])
        for c in ['Chromosome','Start_Position','ref_counts','var_counts','ObsPyCloneCluster','PyClonePhyloCluster']:
            submaf[c] = submaf[c].fillna(0).astype('Int64', errors='ignore')
        #submaf = submaf.astype({'Chromosome':'int','Start_Position':'int','ObsPyCloneCluster':'int','PyClonePhyloCluster':'int'}, errors='ignore')
        submaf = submaf.sort_values(by=['sample_id','Chromosome','Start_Position'], ascending=[True,True,True])
        submaf = submaf.loc[:, submaf.columns.isin(cols)].reindex(columns=cols)
        submaf.to_csv(outdir+'/'+i+'.tsv', sep='\t', index=False, mode='w')


def main():
    ''' '''
    args = GetArgs()
    CombineDfPyclone(args.maf, args.clin, args.cnv, args.out, args.by, args.depth, args.ad, args.af)


if __name__ == '__main__':
    main()
