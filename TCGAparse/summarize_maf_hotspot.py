#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2021/12/31
# @Author  : zhangsiwen
# @contact : zhang.siwen@puruijizhun.com
# @LastModified: 2022/01/05
# @ChangeLog
#     20211231, first version
#     20220105, finished first version


import pandas as pd
import argparse


def GetArgs():
    parser = argparse.ArgumentParser(description='Summarize TCGA tumor hotspot according to maf.', formatter_class=argparse.RawTextHelpFormatter)
    parser._action_groups.pop()
    required = parser.add_argument_group('Required arguments')
    required.add_argument('--maf', help='TCGA maf file, Add.CancerType.maf', action='store', dest='maf', required=True)
    required.add_argument('--info', help='TCGA sample barcode and disease subtype file, Tumor_subtypes_for_PanCanPathways_9125.txt', action='store', dest='info', required=True)
    required.add_argument('--tumor', help='tumor type to be summarized, LUAD', action='store', dest='tumor', required=True)
    required.add_argument('--out', help='output prefix', action='store', dest='out', required=True)
    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument('--by', help='summarize hotspot by "tumor" or "subtype"', choices=['tumor','subtype'], action='store', dest='by', default='subtype', required=False)
    optional.add_argument('--variant', help="Variant_Type type to be reported, separate by comma ',', use 'ALL' for NOT FILTER any type \ndefault 'SNP', choose from [ALL,DEL,INS,SNP,ONP,TNP]", choices=['ALL','DEL','INS','SNP','ONP','TNP'], action='store', dest='variant', default='SNP', required=False)
    optional.add_argument('--func', help="Variant_Classification type to be reported, separate by comma ',', use 'ALL' for NOT FILTER any classification \ndefault 'Frame_Shift_Del,Frame_Shift_Ins,In_Frame_Del,In_Frame_Ins,Intron,Missense_Mutation,Nonsense_Mutation,Nonstop_Mutation,RNA,Silent,Splice_Site,Translation_Start_Site' \nchoose from [ALL,3'Flank,3'UTR,5'Flank,5'UTR,Frame_Shift_Del,Frame_Shift_Ins,In_Frame_Del,In_Frame_Ins,Intron,Missense_Mutation,Nonsense_Mutation,Nonstop_Mutation,RNA,Silent,Splice_Site,Translation_Start_Site]", choices=['ALL','3\'Flank','3\'UTR','5\'Flank','5\'UTR','Frame_Shift_Del','Frame_Shift_Ins','In_Frame_Del','In_Frame_Ins','Intron','Missense_Mutation','Nonsense_Mutation','Nonstop_Mutation','RNA','Silent','Splice_Site','Translation_Start_Site'], action='store', dest='func', default="Frame_Shift_Del,Frame_Shift_Ins,In_Frame_Del,In_Frame_Ins,Intron,Missense_Mutation,Nonsense_Mutation,Nonstop_Mutation,RNA,Silent,Splice_Site,Translation_Start_Site", required=False)
    args = parser.parse_args()
    return args

'''
maffile = 'TCGA.LUAD.maf'
infofile = 'Tumor_subtypes_for_PanCanPathways_9125.txt'
tumor = 'LUAD'
by = 'subtype'
outfile = 'TCGA.LUAD.maf'
func = 'Frame_Shift_Del,Frame_Shift_Ins,In_Frame_Del,In_Frame_Ins,Intron,Missense_Mutation,Nonsense_Mutation,Nonstop_Mutation,RNA,Silent,Splice_Site,Translation_Start_Site'
variant = 'SNP'
info, total_sample = SampleBarcode(infofile, tumor)
maf = FilterMaf(maffile, tumor, func, variant, outfile, total_sample)
'''

def SampleBarcode(infofile, tumor):
    '''
    @ func: read sample barcode and disease file and return info and total sample number
    '''
    info = pd.read_table(infofile, low_memory=False, comment='#')
    info = info[info['DISEASE']==tumor]
    total_sample = info['PATIENT_BARCODE'].nunique()
    print('Total number of sample is %s for %s.' % (total_sample, tumor))
    return info, total_sample


def CheckEmpty(df, condition):
    '''
    @ func: check if dataframe is empty
    '''
    # if not empty, return mutation records number and sample number
    if not df.empty:
        mutation = len(df.index)
        sample = len(set(df['Tumor_Sample_Barcode']))
        print('Filtering %s, %s records of %s sample left.' % (condition, mutation, sample))
        return mutation, sample
    # if empty
    else:
        print('Empty record when filtering %s ' % (condition))
        print('Ending process, please try again with other filter conditions.')
        exit()
    

def FilterMaf(maffile, tumor, func, variant, outfile, total_sample):
    '''
    @ func: filter maf by tumor type, mutant function, variant type, export number summary file
    '''
    # original maf
    maf = pd.read_table(maffile, low_memory=False, comment='#')
    totalMutation, totalMutationSample = CheckEmpty(maf, 'none')
    # filter tumor type
    maf = maf[maf['DISEASE']==tumor]
    containMutation, containMutationSample = CheckEmpty(maf, tumor)
    # filter variant 
    variant = [i.strip() for i in variant.strip().split(',')]
    maf = maf[maf['Variant_Type'].isin(variant)]
    filterVariant, filterVariantSample = CheckEmpty(maf, variant)
    # filter func
    func = [i.strip() for i in func.strip().split(',')]
    maf = maf[maf['Variant_Classification'].isin(func)]
    filterFunc, filterFuncSample = CheckEmpty(maf, func)
    # output record number
    out = open('.'.join([outfile,tumor,'MutationSummary.tsv']), 'w')
    out.write('Type\tMutationNumber\tSampleNumber\n')
    out.write('TotalTumor\t%s\t%s\n' % (totalMutation, total_sample))
    out.write('FilterMutation\t%s\t%s\n' % (containMutation, containMutationSample))
    out.write('FilterVariant\t%s\t%s\n' % (filterVariant, filterVariantSample))
    out.write('FilterFunction\t%s\t%s\n' % (filterFunc, filterFuncSample))
    out.close()
    return maf


def ProcessMaf(submaf, subtype, out):
    '''
    @ func: summarize occurrences, calculate covered sample numbers, for each sites and each genes in submaf
    @ submaf: maf in dataframe format
    @ subtype: string, used in output file name
    '''
    # cols for extracting columns for site comparing
    cols = ['Hugo_Symbol','Chromosome','Start_Position','End_Position','Variant_Classification','Variant_Type','Reference_Allele','HGVSc','HGVCp','Transcript','Existing_variation','RefSeq']
    #################
    # calculate sample number for each site
    #################
    #sites = submaf.pivot_table(columns=cols, aggfunc='size').sort_values(ascending=False)
    sites = submaf.groupby(cols[1:6],as_index=False).size().rename(columns={'size':'Occurrence'}).sort_values(by='Occurrence', ascending=False).reset_index().drop('index', axis=1)
    # how many samples each gene covered, used in output dataframe
    sites_sampleid = []
    sites_sampleid_numbers = []
    sites_sampleid_numbers_cumulate = []
    for s in sites.index:
        # extract sites by index and compare with submaf
        s_sampleid = set(pd.merge(submaf, sites.iloc[[s]], on=cols[1:6], how='inner')['Tumor_Sample_Barcode'])
        sites_sampleid.append(s_sampleid)
        sites_sampleid_numbers.append(len(s_sampleid))
        if sites_sampleid_numbers_cumulate == []:
            sites_sampleid_numbers_cumulate.append(len(s_sampleid))
        else:
            sites_sampleid_numbers_cumulate.append(len(set([j for i in sites_sampleid for j in i])))
    # column-bind DataFrames related to sites into new DataFrame
    sites_matrix = pd.concat([sites, pd.DataFrame({'Covered_Sample':sites_sampleid, 'Covered_Sample_Number':sites_sampleid_numbers, 'Cumulate_Sample_Number':sites_sampleid_numbers_cumulate})], axis=1)
    sites_matrix.to_csv('.'.join([out,subtype,'HotSites.tsv']), sep='\t', index=False)
    #################        
    # calculate sample number for each gene
    #################
    #genes = submaf.pivot_table(columns='Hugo_Symbol', aggfunc='size').sort_values(ascending=False)
    genes = submaf.groupby('Hugo_Symbol',as_index=False).size().rename(columns={'size':'Occurrence'}).sort_values(by='Occurrence', ascending=False).reset_index().drop('index', axis=1)
    # how many samples each gene covered, used in output dataframe
    genes_sampleid = []
    genes_sampleid_numbers = []
    genes_sampleid_numbers_cumulate = []
    for g in genes['Hugo_Symbol']:
        # extract genes by symbol and compare with submaf
        g_samples = set(pd.merge(submaf, genes[genes['Hugo_Symbol']==g], on='Hugo_Symbol', how='inner')['Tumor_Sample_Barcode'])
        genes_sampleid.append(g_samples)
        genes_sampleid_numbers.append(len(g_samples))
        if genes_sampleid_numbers_cumulate == []:
            genes_sampleid_numbers_cumulate.append(len(g_samples))
        else:
            # all sampleid dedup number
            genes_sampleid_numbers_cumulate.append(len(set([j for i in genes_sampleid for j in i])))
    # column-bind DataFrames related to genes into new DataFrame
    genes_matrix = pd.concat([genes, pd.DataFrame({'Covered_Sample':genes_sampleid, 'Covered_Sample_Number':genes_sampleid_numbers, 'Cumulate_Sample_Number':genes_sampleid_numbers_cumulate})], axis=1)
    genes_matrix.to_csv('.'.join([out,subtype,'HotGenes.tsv']), sep='\t', index=False)


def SummarizeMaf(maf, tumor, by, out):
    '''
    @ func: choose to summarize by tumor or subtype
    '''
    # check if subtype applicable
    subtype = list(set(maf['SUBTYPE']))
    nosubtype = subtype == ['Not_Applicable']
    # summarize by subtype
    if by == 'subtype':
        # check if no subtype exists
        if nosubtype:
            print('Subtype not applicable, continue summarizing by tumor type.')
            ProcessMaf(maf, tumor, out)
            print("Finished.")
        # for each subtype
        else:
            output = open('.'.join([out,tumor,'MutationSummary.tsv']), 'a')
            for st in subtype:
                print('Summarizing by %s ...' % (st))
                submaf = maf[maf['SUBTYPE']==st]
                mutation_number, sample_number = CheckEmpty(submaf, st)
                output.write('%s\t%s\t%s\n' % (st, mutation_number, sample_number))
                ProcessMaf(submaf, tumor+"_"+st, out)
                print("Finished.")
            output.close()
    # summarize by tumor
    else:
        print('Summarizing by %s ...' % (tumor))
        ProcessMaf(maf, tumor, out)
        print("Finished.")


def main():
    ''' '''
    args = GetArgs()
    if args.func == 'ALL':
        args.func = "3'Flank,3'UTR,5'Flank,5'UTR,Frame_Shift_Del,Frame_Shift_Ins,In_Frame_Del,In_Frame_Ins,Intron,Missense_Mutation,Nonsense_Mutation,Nonstop_Mutation,RNA,Silent,Splice_Site,Translation_Start_Site"
    if args.variant == 'ALL':
        args.variant = 'DEL,INS,SNP,ONP,TNP'
    # get barcode matrix and total_sample number
    info, total_sample = SampleBarcode(args.info, args.tumor)
    # filter maf by parameters
    maf = FilterMaf(args.maf, args.tumor, args.func, args.variant, args.out, total_sample)
    # summarize maf by hot sites and genes
    SummarizeMaf(maf, args.tumor, args.by, args.out)


if __name__ == '__main__':
    main()
