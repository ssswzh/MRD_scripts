#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2021/12/31
# @Author  : zhangsiwen
# @contact : zhang.siwen@puruijizhun.com
# @LastModified: 2022/01/14
# @ChangeLog
#     20211231, first version
#     20220105, finished first version
#     20220111, add CompareMaf(), modified ProcessMaf(), --base, --mode, --site
#     20220114, add AddColumnPercent(), RenameColumns()


import os
import sys
import pandas as pd
import argparse
import re


def GetArgs():
    parser = argparse.ArgumentParser(description='Summarize tumor hotspot according to maf, compare hotspot according to baseline.', formatter_class=argparse.RawTextHelpFormatter)
    parser._action_groups.pop()
    required = parser.add_argument_group('Required arguments')
    required.add_argument('--maf', help='maf file\nif choose to only compare to baseline (--mode==Compare --base), specify summarized file XX.HotSites.tsv, otherwise specify regular maf file', action='store', dest='maf', required=True)
    required.add_argument('--tumor', help='tumor type to be summarized, LUAD', action='store', dest='tumor', required=True)
    required.add_argument('--out', help='output prefix', action='store', dest='out', required=True)
    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument('--info', help='sample barcode and disease subtype file\nif choose --mode==Summary or --mode==Both, specify info file', action='store', dest='info')
    optional.add_argument('--by', help='summarize hotspot by "Tumor" or "Subtype"', choices=['Tumor','Subtype'], action='store', dest='by', default='subtype', required=False)
    optional.add_argument('--site', help='summarize hotspot by sites in column "SITE", do not specify if you do not have this column', choices=['ALL','Primary','Metastasis','Recurrence','Not_Applicable'], action='store', dest='site', default='ALL', required=False)
    optional.add_argument('--variant', help="Variant_Type type to be reported, separate by comma ',', use 'ALL' for NOT FILTER any type \ndefault 'ALL', choose from [ALL,DEL,INS,SNP,DNP,TNP,ONP]", choices=['ALL','DEL','INS','SNP','DNP','TNP','ONP'], action='store', dest='variant', default='ALL', required=False)
    optional.add_argument('--func', help="Variant_Classification type to be reported, separate by comma ',', use 'ALL' for NOT FILTER any classification \ndefault 'Frame_Shift_Del,Frame_Shift_Ins,In_Frame_Del,In_Frame_Ins,Intron,Missense_Mutation,Nonsense_Mutation,Nonstop_Mutation,RNA,Silent,Splice_Site,Translation_Start_Site' \nchoose from [ALL,3'Flank,3'UTR,5'Flank,5'UTR,Frame_Shift_Del,Frame_Shift_Ins,In_Frame_Del,In_Frame_Ins,Intron,Missense_Mutation,Nonsense_Mutation,Nonstop_Mutation,RNA,Silent,Splice_Site,Translation_Start_Site]", choices=['ALL','3\'Flank','3\'UTR','5\'Flank','5\'UTR','Frame_Shift_Del','Frame_Shift_Ins','In_Frame_Del','In_Frame_Ins','Intron','Missense_Mutation','Nonsense_Mutation','Nonstop_Mutation','RNA','Silent','Splice_Site','Translation_Start_Site'], action='store', dest='func', default="Frame_Shift_Del,Frame_Shift_Ins,In_Frame_Del,In_Frame_Ins,Intron,Missense_Mutation,Nonsense_Mutation,Nonstop_Mutation,RNA,Silent,Splice_Site,Translation_Start_Site", required=False)
    optional.add_argument('--base', help='mutation site baseline maf, if specified, compare input maf to baseline maf\nif --by==subtype, accept key:value pair for baseline maf of each subtype OR one baseline maf for comparing, separate key:value pair by comma ","\nFor types do not have baseline file, use "Other:XXX.tsv" for comparison.\nif --by==tumor, only accept one baseline maf file', action='store', dest='base', required=False)
    optional.add_argument('--mode', help='Summary OR Compare OR Both', choices=['Summary','Compare','Both'], action='store', default='Summary', dest='mode', required=False)
    # usage examples
    usage = '''Usage:
    Summary:
        %(prog)s --maf Add.CancerType.maf --info Tumor_subtypes_for_PanCanPathways_9125.txt --tumor BRCA --out BRCA_hotspot/TCGA --by subtype --variant ALL > log
    Compare:
        %(prog)s --maf TCGA.BRCA.HotSites.tsv --tumor BRCA --out output --by subtype --site Primary --base TCGA.BRCA.HotSites.tsv --mode Compare > log
    Summary and Compare:
        %(prog)s --maf BRCA/brca_ink4_msk_2021/data_mutation.add_clinical.txt --info BRCA/brca_ink4_msk_2021/data_clinical_sample.info.txt --tumor BRCA --out output --by subtype --site Primary --variant ALL --base TCGA.BRCA.HotSites.tsv --mode Both > log
        %(prog)s --maf data_mutation.add_clinical.txt --info data_clinical_sample.info.txt --tumor BRCA --out output --by subtype --site Primary --variant ALL --base Her2:BRCA_Her2.HotSites.tsv,LumA:BRCA_LumA.HotSites.tsv,Other:BRCA.HotSites.tsv --mode Both > log
    '''
    parser.epilog = usage
    args = parser.parse_args()
    return args

'''

maffile = 'BRCA/brca_msk_2018/data_mutation.add_clinical.txt'
infofile = 'BRCA/brca_msk_2018/data_clinical_sample.info.txt'
outfile = "test"

tumor = 'BRCA'
by = 'subtype'
func = 'Frame_Shift_Del,Frame_Shift_Ins,In_Frame_Del,In_Frame_Ins,Intron,Missense_Mutation,Nonsense_Mutation,Nonstop_Mutation,RNA,Silent,Splice_Site,Translation_Start_Site'
variant = 'SNP,DEL,INS,ONP,TNP'
site = 'Primary'

# maffile = '/mnt/ddngs/zhangsw/project/HotDatasets/TCGA_data/BRCA_hotspot/TCGA.BRCA_Basal.HotSites.tsv'
basefile = '/mnt/ddngs/zhangsw/project/HotDatasets/TCGA_data/tumor_BRCA_hotspot/TCGA.BRCA.HotSites.tsv'
info, total_sample = SampleBarcode(infofile, tumor)
maf = FilterMaf(maffile, tumor, func, variant, outfile, total_sample, site)

'''

def SampleBarcode(infofile, tumor):
    '''
    @ func: read sample barcode and disease file and return info and total sample number
    @ invoked by main()
    '''
    info = pd.read_table(infofile, low_memory=False, comment='#')
    info = info[info['DISEASE']==tumor]
    total_sample = info['PATIENT_BARCODE'].nunique()
    print('Total number of sample is %s for %s.' % (total_sample, tumor))
    return info, total_sample


def CheckEmpty(df, condition):
    '''
    @ func: check if dataframe is empty
    @ invoked by FilterMaf()
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
        sys.exit('Ending process, please try again with other filter conditions.')
    

def FilterMaf(maffile, tumor, func, variant, outfile, total_sample, site):
    '''
    @ func: filter maf by tumor type, mutant function, variant type, export number summary file
    @ invoked by main()
    '''
    # original maf
    maf = pd.read_table(maffile, low_memory=False, comment='#')
    totalMutation, totalMutationSample = CheckEmpty(maf, 'none')
    # filter tumor type
    maf = maf[maf['DISEASE']==tumor]
    containMutation, containMutationSample = CheckEmpty(maf, tumor)
    # if choose to filter site (eg. Primary, Metastasis, Reccurence, Not_Applicable)
    if site:
        site = [i.strip() for i in site.strip().split(',')]
        maf = maf[maf['SITE'].isin(site)]
        filterSite, filterSiteSample = CheckEmpty(maf, site)
    # filter variant 
    if variant:
        variant = [i.strip() for i in variant.strip().split(',')]
        maf = maf[maf['Variant_Type'].isin(variant)]
        filterVariant, filterVariantSample = CheckEmpty(maf, variant)
    # filter func
    if func:
        func = [i.strip() for i in func.strip().split(',')]
        maf = maf[maf['Variant_Classification'].isin(func)]
        filterFunc, filterFuncSample = CheckEmpty(maf, func)
    # output record number
    out = open('.'.join([outfile,tumor,'MutationSummary.tsv']), 'w')
    out.write('Type\tMutationNumber\tSampleNumber\n')
    out.write('TotalTumor\t%s\t%s\n' % (totalMutation, total_sample))
    out.write('FilterMutation\t%s\t%s\n' % (containMutation, containMutationSample))
    if site:
        out.write('FilterSite\t%s\t%s\n' % (filterSite, filterSiteSample))
    if variant:
        out.write('FilterVariant\t%s\t%s\n' % (filterVariant, filterVariantSample))
    if func:
        out.write('FilterFunction\t%s\t%s\n' % (filterFunc, filterFuncSample))
    out.close()
    return maf


def AddColumnPercent(df, total, col, name):
    '''
    @ func: add percent column to dataframe
    '''
    # calculate cumulate sample percent for each site
    df[name] = df[col]/total
    return df


def SummarizeMaf(submaf, out, command):
    '''
    @ func: summarize occurrences, calculate covered sample numbers, for each sites and each genes in submaf
    @ invoked by ProcessMaf()
    @ submaf: maf in dataframe format
    '''
    # cols for extracting columns for site comparing
    cols = ['Hugo_Symbol','Chromosome','Start_Position','End_Position','Variant_Classification','Variant_Type','Reference_Allele','HGVSc','HGVCp','Transcript','Existing_variation','RefSeq']
    #################
    # calculate sample number for each site
    #################
    #sites = submaf.pivot_table(columns=cols, aggfunc='size').sort_values(ascending=False)
    sites = submaf.groupby(cols[0:6],as_index=False).size().rename(columns={'size':'Occurrence'}).sort_values(by='Occurrence', ascending=False).reset_index().drop('index', axis=1)
    # how many samples each gene covered, used in output dataframe
    sites_sampleid = []
    sites_sampleid_numbers = []
    sites_sampleid_numbers_cumulate = []
    for s in sites.index:
        # extract sites by index and compare with submaf
        s_sampleid = set(pd.merge(submaf, sites.iloc[[s]], on=cols[0:6], how='inner')['Tumor_Sample_Barcode'])
        sites_sampleid.append(s_sampleid)
        sites_sampleid_numbers.append(len(s_sampleid))
        if sites_sampleid_numbers_cumulate == []:
            sites_sampleid_numbers_cumulate.append(len(s_sampleid))
        else:
            sites_sampleid_numbers_cumulate.append(len(set([j for i in sites_sampleid for j in i])))
    # column-bind DataFrames related to sites into new DataFrame
    sites_matrix = pd.concat([sites, pd.DataFrame({'Covered_Sample':sites_sampleid, 'Covered_Sample_Number':sites_sampleid_numbers, 'Cumulate_Sample_Number':sites_sampleid_numbers_cumulate})], axis=1)
    sites_matrix = AddColumnPercent(sites_matrix, max(sites_matrix['Cumulate_Sample_Number']), 'Covered_Sample_Number', 'Covered_Sample_Number_Percent')
    sites_matrix = AddColumnPercent(sites_matrix, max(sites_matrix['Cumulate_Sample_Number']), 'Cumulate_Sample_Number', 'Cumulate_Sample_Number_Percent')
    out_site = '.'.join([out,'HotSites.tsv'])
    with open(out_site, 'w') as f:
        f.write('# Command: %s\n' % command)    
    sites_matrix.to_csv(out_site, sep='\t', index=False, mode='a')
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
    genes_matrix = AddColumnPercent(genes_matrix, max(genes_matrix['Cumulate_Sample_Number']), 'Covered_Sample_Number', 'Covered_Sample_Number_Percent')
    genes_matrix = AddColumnPercent(genes_matrix, max(genes_matrix['Cumulate_Sample_Number']), 'Cumulate_Sample_Number', 'Cumulate_Sample_Number_Percent')
    out_gene = '.'.join([out,'HotGenes.tsv'])
    with open(out_gene, 'w') as f:
        f.write('# Command: %s\n' % command)
    genes_matrix.to_csv(out_gene, sep='\t', index=False, mode='a')
    return sites_matrix, genes_matrix


def RenameColumns(columns, prefix):
    '''
    @ func: for comparing, rename dataframe column names by adding prefix to specific columns
    @ invoked by CompareMaf()
    '''
    colrename = list(columns)
    pattern = ['Occurrence', 'Covered', 'Cumulate']
    for c in range(len(colrename)):
        if colrename[c].startswith(tuple(pattern)):
            colrename[c] = prefix + colrename[c]
    return colrename


def CompareMaf(maf, basefile, out, command):
    '''
    @ func: compare maf
    @ invoked by ProcessMaf()
    '''
    # read baseline file
    base = pd.read_table(basefile, low_memory=False, comment='#')
    # cols for extracting columns for site comparing
    cols = ['Hugo_Symbol','Chromosome','Start_Position','End_Position','Variant_Classification','Variant_Type','Reference_Allele','HGVSc','HGVCp','Transcript','Existing_variation','RefSeq']
    # calculate cumulate sample percent for each site
    base = AddColumnPercent(base, max(base['Cumulate_Sample_Number']), 'Cumulate_Sample_Number', 'Cumulate_Sample_Number_Percent')
    maf = AddColumnPercent(maf, max(maf['Cumulate_Sample_Number']), 'Cumulate_Sample_Number', 'Cumulate_Sample_Number_Percent')
    # change column names, preparing for dataframe merge
    base.columns = RenameColumns(base.columns, 'Baseline_')
    maf.columns = RenameColumns(maf.columns, 'Dataset_')
    # merge baseline and maf, fill NaN by zero, delete dataset cumulate sample number and percent
    #mergedDf = base.merge(maf, on=cols[1:6], how='outer').sort_values(by=['Baseline_Occurrence','Dataset_Occurrence'], ascending=[False,False]).reset_index().drop('index', axis=1).fillna('')
    mergedDf = base.merge(maf, on=cols[1:6], how='outer').fillna('')
    del mergedDf['Dataset_Cumulate_Sample_Number']
    del mergedDf['Dataset_Cumulate_Sample_Number_Percent']
    # re-count sample coverage for each site in merged dataframe
    merged_sampleid = []
    merged_sampleid_numbers_cumulate = []
    for s in mergedDf.index:
        if mergedDf.iloc[s,]['Dataset_Covered_Sample'] != '':
            # extract sites by index and compare with submaf
            # get sample ids
            covered_sample = mergedDf.iloc[s,]['Dataset_Covered_Sample']
            if type(covered_sample) is set:
                merged_sampleid.append(covered_sample)
            else:
                # sample id is string type in 'set()' display format
                s_sampleid = set([i.strip() for i in re.sub("['{}]","",covered_sample).split(',')])
                merged_sampleid.append(s_sampleid)
            if merged_sampleid_numbers_cumulate == []:
                merged_sampleid_numbers_cumulate.append(len(merged_sampleid[0]))
            else:
                merged_sampleid_numbers_cumulate.append(len(set([j for i in merged_sampleid for j in i])))
        else:
            merged_sampleid.append('')
            if merged_sampleid_numbers_cumulate == []:
                merged_sampleid_numbers_cumulate.append(0)
            else:
                merged_sampleid_numbers_cumulate.append(merged_sampleid_numbers_cumulate[-1])
    # add re-count sample coverage and percent to new DataFrame
    mergedDf = pd.concat([mergedDf, pd.DataFrame({'Dataset_Cumulate_Sample_Number':merged_sampleid_numbers_cumulate})], axis=1)
    mergedDf = AddColumnPercent(mergedDf, max(mergedDf['Dataset_Cumulate_Sample_Number']), 'Dataset_Cumulate_Sample_Number', 'Dataset_Cumulate_Sample_Number_Percent')
    # write out merged dataframe and comment
    outfile = '.'.join([out,'HotSites.tsv'])
    with open(outfile, 'w') as f:
        f.write('# Command: %s\n' % command)
    mergedDf.to_csv(outfile, sep='\t', index=False, mode='a')


def ProcessMaf(maf, tumor, by, out, base, mode, command):
    '''
    @ func: choose to summarize by tumor or subtype
    @ invoked by main()
    '''
    # check how many baseline files are included
    if len(base.strip().split(',')) == 1:
        base = base.strip()
    else: # convert to dict
        base = {i.strip().split(':')[0].strip():i.strip().split(':')[1].strip() for i in base.strip().split(',')}
    # summarize and compare by subtype
    if by == 'Subtype':
        # check if subtype applicable
        subtype = list(set(maf['SUBTYPE']))
        nosubtype = subtype == ['Not_Applicable']
        # check if no subtype exists
        if nosubtype:
            if type(base) is dict:
                sys.exit('Subtype not applicable, please specify only one --base file.')
            if mode == 'Summary':
                print('Subtype not applicable, continue summarizing by tumor type.')
                sites_matrix, genes_matrix = SummarizeMaf(maf, '.'.join([out,'Summary',tumor]), command)
            if mode == 'Compare':
                print('Subtype not applicable, continue comparing by tumor type.')
                CompareMaf(maf, base, '.'.join([out,'Compare',tumor]), command)
            if mode == 'Both':
                print('Subtype not applicable, continue summarizing by tumor type.')
                sites_matrix, genes_matrix = SummarizeMaf(maf, '.'.join([out,'Summary',tumor]), command)
                print('Subtype not applicable, continue comparing by tumor type.')
                CompareMaf(sites_matrix, base, '.'.join([out,'Compare',tumor]), command)
        # for each subtype
        else:
            output = open('.'.join([out,tumor,'MutationSummary.tsv']), 'a')
            for st in subtype:
                submaf = maf[maf['SUBTYPE']==st]
                mutation_number, sample_number = CheckEmpty(submaf, st)
                output.write('%s\t%s\t%s\n' % (st, mutation_number, sample_number))
                if mode == 'Summary':
                    print('Summarizing by %s ...' % (st))
                    sites_matrix, genes_matrix = SummarizeMaf(submaf, '.'.join([out,'Summary',tumor+"_"+st]), command)
                if mode == 'Compare':
                    print('Comparing by %s ...' % (st))
                    # if specify multiple baseline files
                    if len(base) > 1:
                        if st in base:
                            CompareMaf(submaf, base[st], '.'.join([out,'Compare',tumor+"_"+st]), command)
                        elif 'Other' in base:
                            CompareMaf(submaf, base['Other'], '.'.join([out,'Compare',tumor+"_"+st]), command)
                        else:
                            print('Subtype %s do not have corresponding baseline file or "Other:XXXfile" for comparison, skipping...' % st)
                            print('You can use "Other:XXXfile" for comparing subtypes that do not have specific baseline file.')
                    else:
                        CompareMaf(submaf, base, '.'.join([out,'Compare',tumor+"_"+st]), command)
                if mode == 'Both':
                    print('Summarizing by %s ...' % (st))
                    sites_matrix, genes_matrix = SummarizeMaf(submaf, '.'.join([out,'Summary',tumor+"_"+st]), command)
                    print('Comparing by %s ...' % (st))
                    if len(base) > 1:
                        if st in base:
                            CompareMaf(sites_matrix, base[st], '.'.join([out,'Compare',tumor+"_"+st]), command)
                        elif 'Other' in base:
                            CompareMaf(sites_matrix, base['Other'], '.'.join([out,'Compare',tumor+"_"+st]), command)
                        else:
                            print('Subtype %s do not have corresponding baseline file or "Other:XXXfile" for comparison, skipping...' % st)
                            print('You can use "Other:XXXfile" for comparing subtypes that do not have specific baseline file.')
                    else:
                        CompareMaf(sites_matrix, base, '.'.join([out,'Compare',tumor+"_"+st]), command)
            output.close()
    # summarize and compare by tumor
    else:
        if mode == 'Summary':
            print('Summarizing by %s ...' % (tumor))
            sites_matrix, genes_matrix = SummarizeMaf(maf, '.'.join([out,'Summary',tumor]), command)
        if mode == 'Compare':
            print('Comparing by %s ...' % (tumor))
            CompareMaf(maf, base, '.'.join([out,'Compare',tumor]), command)
        if mode == 'Both':
            print('Summarizing by %s ...' % (tumor))
            sites_matrix, genes_matrix = SummarizeMaf(maf, '.'.join([out,'Summary',tumor]), command)
            print('Comparing by %s ...' % (tumor))
            CompareMaf(sites_matrix, base, '.'.join([out,'Compare',tumor]), command)
    print("Finished.")


def main():
    ''' '''
    args = GetArgs()
    # command line arguments
    command = [os.path.basename(__file__)]
    for k, v in dict(args._get_kwargs()).items():
        if v != None:
            command.append('--'+k)
            command.append(v)
    command = ' '.join(command)
    # dependencies
    if (args.mode=='Compare' or args.mode=='Both') and args.base is None:
        sys.exit('If choose --mode==Compare or --mode==Both, please also specify --base file')
    if args.by=='Tumor' and len(args.base.strip().split(','))!=1:
        sys.exit('If choose --by==Tumor, please specify a single --base file')
    if (args.mode=='Summary' or args.mode=='Both') and args.info is None:
        sys.exit('If choose --mode==Summary or --mode==Both, please specify --info file')
    # assign values
    if args.func == 'ALL':
        args.func = ''
    if args.variant == 'ALL':
        args.variant = ''
    if args.site == 'ALL':
        args.site = ''
    # if choose to summary or both, filter maf by parameters
    if args.mode=='Summary' or args.mode=='Both':
        # get barcode matrix and total_sample number
        info, total_sample = SampleBarcode(args.info, args.tumor)
        maf = FilterMaf(args.maf, args.tumor, args.func, args.variant, args.out, total_sample, args.site)
    else: # else (for compare) only read maf (HotSites.tsv)
        maf = pd.read_table(args.maf, low_memory=False, comment='#')
        if 'Covered_Sample' not in maf.columns:
            sys.exit('Error: for comparing-only purpose, maf file should be HotSites.tsv generated by this program in Summary mode.')
    # summarize and compare maf by hot sites and genes
    ProcessMaf(maf, args.tumor, args.by, args.out, args.base, args.mode, command)


if __name__ == '__main__':
    main()



