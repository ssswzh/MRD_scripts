#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2022/03/07
# @Author  : zhangsiwen
# @contact : zhang.siwen
# @LastModified: 2022/08/24
# @ChangeLog
#     20220309, first version
#     20220311, add --hotspot
#     20220315, check if total number larger than mutation number
#     20220609, fix problem that sometimes the total number of selection is less than the given number due to GC content
#     20220630, add --maf and AnnotateSites() to annotate final site selection result
#     20220818, add --af, ignore sites with VAF<5% (which are probabily undetectable)
#     20220824, re-write the whole script, version 2


import sys
import argparse
import pandas as pd
import numpy as np
from Bio.SeqUtils import GC 
from pyfaidx import Fasta
from Bio.Seq import Seq


def GetArgs():
    parser = argparse.ArgumentParser(description='Select sites from pyclone result', formatter_class=argparse.RawTextHelpFormatter)
    parser._action_groups.pop()
    # required arguments
    required = parser.add_argument_group('Required arguments')
    required.add_argument('--tsv', help='pyclone result h5.tsv', action='store', dest='tsv', required=True)
    required.add_argument('--dat', help='pyclone input file', action='store', dest='dat', required=True)
    required.add_argument('--out', help='out file name', action='store', dest='out', required=True)
    # optional arguments
    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument('--ref', help='reference fasta file', default='/mnt/user/errand/pipeline/MC3/ref/Homo_sapiens_assembly19.fasta', action='store', dest='ref', required=False)
    optional.add_argument('--site', help='total site number of this selection, default [30]', default=30, type=int, action='store', dest='site', required=False)
    optional.add_argument('--af', help='filter variant with allele frequency larger than this, default [0.05]', default=0.05, type=float, action='store', dest='af', required=False)
    optional.add_argument('--gc', help='GC content range separated by comma ",", default [30,70]', default='[30,70]', action='store', dest='gc', required=False)
    optional.add_argument('--flank', help='flank sequence length of both sides, default [50]', default=50, type=int, action='store', dest='flank', required=False)
    optional.add_argument('--min', help='minimum number of site for each cluster if its proportion is too low, default [3]', default=3, type=int, action='store', dest='min', required=False)
    optional.add_argument('--hotspot', help='if give hotspot bed file, the mutation located in regions will be chosen regardless of cluster, \ne.g /mnt/ddngs/zhangsw/project/MRD/OV_WES/hotspots/oncogene_whitelist.annotationCombined.filtered.tsv ', default='', action='store', dest='hotspot', required=False)
    optional.add_argument('--maf', help='provide maf file if you want to annotate sites', action='store', dest='maf', required=False)
    optional.add_argument('--prev', help='if give maf, set population prevalence cutoff in databases, \nFormat: column_name_suffix:cutoff, e.g. AF:0.01, default not filter [AF:1]', default='AF:1', action='store', dest='prev', required=False)
    optional.add_argument('--indel', help='if give maf, filter indel with size larger than this, default [2]', default=2, type=int, action='store', dest='indel', required=False)
    optional.add_argument('--func', help='if give maf, the level to filter mutations functions, default [high] \nlow: Frame_Shift_Del,Frame_Shift_Ins,Missense_Mutation,Nonsense_Mutation,Nonstop_Mutation,Silent,Splice_Site,Translation_Start_Site, \nhigh: Missense_Mutation,Nonsense_Mutation,Nonstop_Mutation, \nUse "None" for not filtering.', default='high' ,action='store', dest='func', required=False)
    #optional.add_argument('--func', help="Variant_Classification type in maf file to be reported, separate by comma ',', use 'ALL' for NOT FILTER any classification \ndefault 'Frame_Shift_Del,Frame_Shift_Ins,In_Frame_Del,In_Frame_Ins,Intron,Missense_Mutation,Nonsense_Mutation,Nonstop_Mutation,RNA,Silent,Splice_Site,Translation_Start_Site' \nchoose from [ALL,3'Flank,3'UTR,5'Flank,5'UTR,Frame_Shift_Del,Frame_Shift_Ins,In_Frame_Del,In_Frame_Ins,Intron,Missense_Mutation,Nonsense_Mutation,Nonstop_Mutation,RNA,Silent,Splice_Site,Translation_Start_Site]", choices=['ALL','3\'Flank','3\'UTR','5\'Flank','5\'UTR','Frame_Shift_Del','Frame_Shift_Ins','In_Frame_Del','In_Frame_Ins','Intron','Missense_Mutation','Nonsense_Mutation','Nonstop_Mutation','RNA','Silent','Splice_Site','Translation_Start_Site'], action='store', dest='func', default="Frame_Shift_Del,Frame_Shift_Ins,In_Frame_Del,In_Frame_Ins,Intron,Missense_Mutation,Nonsense_Mutation,Nonstop_Mutation,RNA,Silent,Splice_Site,Translation_Start_Site", required=False)
    usage = '''Usage:
    %(prog)s --tsv pyclone.h5.tsv --dat pyclone.dat.tsv --out out.tsv
    %(prog)s --tsv pyclone.h5.tsv --dat pyclone.dat.tsv --out out.tsv --ref ref.fasta --site 30 --af 0.05 --gc 30,70 --flank 50 --min 3 --hotspot /mnt/ddngs/zhangsw/project/MRD/OV_WES/hotspots/TCGA_merged.hotspot.tsv --maf mutation.maf --prev AF:0.01 --indel 2
'''
    parser.epilog = usage
    args = parser.parse_args()
    return args


'''
tsv = 'P93.new_pyclone.h5.tsv'
dat = 'P93.new_pyclone.tsv'
out = 'test.tsv'
reference = '/mnt/user/errand/pipeline/MC3/ref/Homo_sapiens_assembly19.fasta'
number = 50
gc = '30,70'
flank = 150
minsite = 3
hotspot = '/mnt/ddngs/zhangsw/project/MRD/OV_WES/hotspots/TCGA_merged.hotspot.tsv'
maf = ''
prev = 'AF:0.01'
indel = 2
'''

# some fixed parameters
#Variant_Classification = ["3'UTR","5'UTR","Frame_Shift_Del","Frame_Shift_Ins","IGR","In_Frame_Del","In_Frame_Ins","Intron","Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation","RNA","Silent","Splice_Region","Splice_Site","Translation_Start_Site"]
Variant_Classification_low = ["Frame_Shift_Del","Frame_Shift_Ins","Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation","Silent","Splice_Site","Translation_Start_Site"]
Variant_Classification_high = ["Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation"]
#Consequence = ['3_prime_UTR_variant','5_prime_UTR_variant','coding_sequence_variant','downstream_gene_variant','frameshift_variant','inframe_deletion','inframe_insertion','intergenic_variant','intron_variant','missense_variant','non_coding_transcript_exon_variant','protein_altering_variant','splice_acceptor_variant','splice_donor_variant','splice_region_variant','start_lost','stop_gained','stop_lost','stop_retained_variant','synonymous_variant','upstream_gene_variant','non_coding_transcript_variant','start_retained_variant']
Consequence_low = ['coding_sequence_variant','frameshift_variant','missense_variant','protein_altering_variant','splice_acceptor_variant','splice_donor_variant','start_lost','stop_gained','stop_lost','stop_retained_variant','synonymous_variant','start_retained_variant']
Consequence_high = ['missense_variant','protein_altering_variant','start_lost','stop_gained','stop_lost']


def GetClusterFraction(df):
    '''
    @ func: get cluster site number and fraction
    @ result: {0: [3, 0.0152], 1: [8, 0.0406], 2: [5, 0.0254], 3: [1, 0.0051], 4: [23, 0.1168], 5: [157, 0.797]}
    '''
    cluster_dict = {}
    total_sites = df.shape[0]
    for i in set(df['cluster_id']):
        num = list(df['cluster_id']).count(i)
        cluster_dict[i] = [num, round(float(num/total_sites),4)]
    # sort dict by value (site number in this case)
    cluster_dict = {r:cluster_dict[r] for r in sorted(cluster_dict, key=cluster_dict.get, reverse=False)}
    return cluster_dict


def RetrieveHotspot(hotspot, sitesDf):
    '''
    @ func: retrieve hotspot sites from provided bed file
    @ result: list of index in sitesDf
    '''
    hotspotDf = pd.read_csv(hotspot, low_memory=False, sep='\t', comment='#', header=0)
    intersectDf = sitesDf.merge(hotspotDf, on=['chrom','pos'], how='inner').fillna('')
    # record intersect index and resultDf
    hotspot_index = []
    for i in intersectDf['mutation_id']:
        hotspot_index.append(sitesDf[sitesDf['mutation_id']==i].index.tolist()[0])
    return hotspot_index


def ExpectClusterSize(clusters, keys, minsite, totalnumber, expect_cluster):
    '''
    @ func: get cluster sizes according to fraction, total number and minimum sites
    @ structure: clusters has the same structure as result of GetClusterFraction()
    @ result: expect_cluster = {3: 1, 0: 3, 1: 3, 2: 3, 4: 6, 5: 34}
    '''
    restnumber = totalnumber
    for k in keys:
        n, p = clusters[k]
        if k not in expect_cluster.keys():
            expect_cluster[k] = 0
        # expected sampling site number is less than minimum site
        if p*totalnumber < minsite:
            # original site number is less than minimum site
            if n < minsite:
                expect_cluster[k] += int(min(restnumber, n))
                restnumber -= int(min(restnumber, n))
            else:
                expect_cluster[k] += int(min(restnumber, minsite))
                restnumber -= int(min(restnumber, minsite))
        # expected sampling site number is larger than minimum site
        else:
            expect_cluster[k] += int(min(restnumber, p*totalnumber))
            restnumber -= int(min(restnumber, p*totalnumber))
    return expect_cluster, restnumber


def MutationFunctionFilter(mafDf, functions):
    ''' 20220824
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


def AnnotateSites(siteDf, maffile, prev_cutoff='AF:1', func='high'):
    ''' 20220630
    @ func: merge sites and maf by mutation_id, return merged dataframe
    '''
    mafDf = pd.read_csv(maffile, low_memory=False, sep='\t', comment='#', header=0)
    if 'CHROM' in mafDf.columns and 'POS' in mafDf.columns:
        mafDf['mutation_id'] = mafDf[['CHROM','POS','Reference_Allele','Tumor_Seq_Allele2']].astype(str).agg(':'.join, axis=1)
    elif 'Chromosome' in mafDf.columns and 'Start_Position' in mafDf.columns:
        mafDf['mutation_id'] = mafDf[['Chromosome','Start_Position','Reference_Allele','Tumor_Seq_Allele2']].astype(str).agg(':'.join, axis=1)
    else:
        sys.exit("Please provide standard maf format which has columns: ['Chromosome','Start_Position'].")
    siteDf = siteDf.merge(mafDf, on='mutation_id', how='inner').fillna('')
    
    # -- filter population prevalence in database, 20220818 --
    if prev_cutoff != '':
        prev_str, prev_af = prev_cutoff.strip().split(':')
        database_cols = [i for i in list(siteDf.columns) if i.endswith(prev_str)]
        for col in database_cols:
            siteDf = siteDf[(siteDf[col]<=float(prev_af)) | (siteDf[col].isnull())]
    
    # -- filter mutation functions --
    if func == 'high':
        functions = [Variant_Classification_high, Consequence_high]
    elif func == 'low':
        functions = [Variant_Classification_low, Consequence_low]
    elif func == 'None':
        return mafDf
    mafDf = MutationFunctionFilter(mafDf, functions)
    
    return siteDf


def SiteSelection(tsv, dat, out, reference, number=30, gc='30,70', af=0.05, flank=50, minsite=3, hotspot='', maffile='', prev_cutoff='AF:1', indel=2, func='high'):
    '''
    @ func: select sites according to cluster fraction, site number, vaf, gc content
    '''
    # -- load reference data --
    # load reference fasta
    fasta = Fasta(reference, rebuild=False)
    
    # -- read files and process --
    pycloneDf = pd.read_csv(tsv, low_memory=False, sep='\t', comment='#', header=0)
    dataDf = pd.read_csv(dat, low_memory=False, sep='\t', comment='#', header=0)
    sitesDf = pycloneDf.merge(dataDf, on='mutation_id', how='outer').fillna('')
    # split 'mutation_id' into four columns
    sitesDf[['chrom','pos','ref','alt']] = sitesDf['mutation_id'].str.split(':', expand=True)
    sitesDf['pos'] = sitesDf['pos'].astype(int)

    # -- flank sequence and gc content --
    # for each site get flank seq and gc content
    for idx in sitesDf.index:
        mutation = sitesDf['mutation_id'][idx]
        chrom, pos, ref, alt = mutation.split(':')
        pos = int(pos)
        # extract sequence from reference and calculate GC content
        flankseq = fasta[chrom][pos-flank-1 : pos+flank].seq
        sitesDf.loc[idx, 'flank_sequence'] = flankseq
        sitesDf.loc[idx, 'gc_content'] = round(GC(flankseq),2)
    
    # -- annotate sites --
    if maffile != '':
        mafflag = True
        sitesDf = AnnotateSites(sitesDf, maffile, prev_cutoff, func)
        # -- indel length --
        for idx in sitesDf.index:
            sitesDf.loc[idx, 'indel_len'] = max(len(sitesDf.loc[idx, 'Reference_Allele']), len(sitesDf.loc[idx, 'Tumor_Seq_Allele2']))
    
    # -- filter sites --
    lowerGC, higherGC = [float(i) for i in gc.strip('[]').split(',')]
    if mafflag:
        filteredDf = sitesDf[(sitesDf['gc_content'] >= lowerGC) & (sitesDf['gc_content'] <= higherGC) & (sitesDf['vaf'] >= af) & (sitesDf['indel_len'] <= indel)]
        restDf = sitesDf[(sitesDf['gc_content'] >= lowerGC) | (sitesDf['gc_content'] <= higherGC) | (sitesDf['vaf'] >= af) | (sitesDf['indel_len'] <= indel)]
    else:
        filteredDf = sitesDf[(sitesDf['gc_content'] >= lowerGC) & (sitesDf['gc_content'] <= higherGC) & (sitesDf['vaf'] >= af)]
        restDf = sitesDf[(sitesDf['gc_content'] >= lowerGC) | (sitesDf['gc_content'] <= higherGC) | (sitesDf['vaf'] >= af)]
    
    # -- output sites after first filtering -- 
    # check if total number larger than mutation number, if true, fill sites to meet the number, 20220609
    mutation_number = filteredDf.shape[0]
    if number > mutation_number:
        # add weight to vaf, gc, indel
        if mafflag:
            restDf['score'] = np.exp(restDf['vaf']*10) + np.log(1/abs(restDf['gc_content']-50)) + np.log(1/restDf['indel_len'])
            del resultDf['indel_len']
        else:
            restDf['score'] = np.exp(restDf['vaf']*10) + np.log(1/abs(restDf['gc_content']-50))
        restDf = restDf.sort_values(by='score', ascending=False).iloc[0:(number-mutation_number), ]
        del resultDf['score']
        resultDf = pd.concat([filteredDf, restDf], axis=0).sort_index()
        resultDf.to_csv(out, sep='\t', index=False, mode='w')
        return None
    elif number == mutation_number:
        if mafflag:
            del resultDf['indel_len']
        filteredDf.to_csv(out, sep='\t', index=False, mode='w')
        return None

    # -- keep hotspot sites from sitesDf --
    if hotspot != '':
        hotspot_index = RetrieveHotspot(hotspot, sitesDf)
        number = number - len(hotspot_index)
        filteredDf = filteredDf.loc[[i for i in filteredDf.index.tolist() if i not in hotspot_index]]

    # -- calculate expect site number for each cluster --
    # cluster id dict with site number and fraction
    filtered_cluster = GetClusterFraction(filteredDf)
    expect_cluster = {}
    restnumber = number
    # get clusters, in case the total number of selected sites not enough
    while True:
        expect_cluster, restnumber = ExpectClusterSize(filtered_cluster, filtered_cluster.keys(), minsite, restnumber, expect_cluster)
        if restnumber == 0:
            break

    # -- get sites --
    idx = []
    for cluster,num in expect_cluster.items():
        subdf = filteredDf[filteredDf['cluster_id'] == cluster].sort_values(by='vaf', ascending=False)
        idx.append(list(subdf.iloc[0:num].index))
    # flatten idx list
    if hotspot != '':
        idx = [j for i in idx for j in i] + hotspot_index
    else:
        idx = [j for i in idx for j in i]
    idx.sort()
    
    # -- write output --
    resultDf = sitesDf.loc[idx].sort_values(by='vaf', ascending=False)
    del resultDf['indel_len']
    resultDf.to_csv(out, sep='\t', index=False, mode='w')
    return None


def main():
    ''' '''
    args = GetArgs()
    if args.maf == '':
        print("Warning: Not providing MAF file, ignore --prev, --indel, --func.")
    SiteSelection(args.tsv, args.dat, args.ref, args.out, args.site, args.gc, args.af, args.flank, args.min, args.hotspot, args.maf, args.prev, args.indel, args.func)


if __name__ == '__main__':
    main()
