#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2022/03/07
# @Author  : zhangsiwen
# @contact : zhang.siwen@puruijizhun.com
# @LastModified: 2022/06/30
# @ChangeLog
#     20220309, first version
#     20220311, add --hotspot
#     20220315, check if total number larger than mutation number
#     20220609, fix problem that sometimes the total number of selection is less than the given number due to GC content
#     20220630, add --maf and AnnotateSites() to annotate final site selection result
#     20220818, add --af, ignore sites with VAF<5% (which are probabily undetectable)


import sys
import argparse
import pandas as pd
from Bio.SeqUtils import GC 
from pyfaidx import Fasta
from Bio.Seq import Seq


def GetArgs():
    parser = argparse.ArgumentParser(description='Select sites from pyclone result', formatter_class=argparse.RawTextHelpFormatter)
    parser._action_groups.pop()
    required = parser.add_argument_group('Required arguments')
    required.add_argument('--tsv', help='pyclone result h5.tsv', action='store', dest='tsv', required=True)
    required.add_argument('--dat', help='pyclone input file', action='store', dest='dat', required=True)
    required.add_argument('--ref', help='reference fasta file', action='store', dest='ref', required=True)
    required.add_argument('--out', help='out file name', action='store', dest='out', required=True)
    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument('--site', help='total site number of this selection, default [30]', default=30, type=int, action='store', dest='site', required=False)
    optional.add_argument('--gc', help='GC content range separated by comma ",", default [30,70]', default='[30,70]', action='store', dest='gc', required=False)
    optional.add_argument('--af', help='filter variant with allele frequency larger than this, default [0.05]', default=0.05, type=float, action='store', dest='af', required=False)
    optional.add_argument('--flank', help='flank sequence length of both sides, default [100]', default=100, type=int, action='store', dest='flank', required=False)
    optional.add_argument('--min', help='minimum number of site for each cluster if its proportion is too low, default [3]', default=3, type=int, action='store', dest='min', required=False)
    optional.add_argument('--hotspot', help='if give hotspot bed file, the mutation located in regions will be chosen regardless of cluster', default='', action='store', dest='hotspot', required=False)
    optional.add_argument('--maf', help='provide maf file if you want to annotate sites', action='store', dest='maf', required=False)
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
'''


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
        # expected sampling site number is larger than minimum site, 
        else:
            expect_cluster[k] += int(min(restnumber, p*totalnumber))
            restnumber -= int(min(restnumber, p*totalnumber))
    return expect_cluster, restnumber


def AnnotateSites(siteDf, maffile):
    # 20220630
    '''
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
    return siteDf


def SiteSelection(tsv, dat, reference, out, number = 50, gc = '30,70', af = 0.05, flank = 150, minsite = 3, hotspot = '', maffile = ''):
    '''
    @ func: select sites according to cluster fraction, site number, vaf, gc content
    '''
    # -- prepare data --
    # load reference fasta
    fasta = Fasta(reference, rebuild=False)
    # GC content
    lowerGC, higherGC = [float(i) for i in gc.strip('[]').split(',')]
    # read pyclone result and input file
    pycloneDf = pd.read_csv(tsv, low_memory=False, sep='\t', comment='#', header=0)
    dataDf = pd.read_csv(dat, low_memory=False, sep='\t', comment='#', header=0)
    sitesDf = pycloneDf.merge(dataDf, on='mutation_id', how='outer').fillna('')
    # split 'mutation_id' into four columns
    sitesDf[['chrom','pos','ref','alt']] = sitesDf['mutation_id'].str.split(':', expand=True)
    sitesDf['pos'] = sitesDf['pos'].astype(int)
    sitesDf = sitesDf[sitesDf['vaf'] >= af]

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
    # filter sites by gc content and vaf
    filteredDf = sitesDf[(sitesDf['gc_content'] >= lowerGC) & (sitesDf['gc_content'] <= higherGC) ]
    # check if total number larger than mutation number, if true, fill sites to meet the number, 20220609
    mutation_number = filteredDf.shape[0]
    if number > mutation_number:
        restDf = sitesDf[(sitesDf['gc_content'] < lowerGC) | (sitesDf['gc_content'] > higherGC)]
        restDf['dif'] = abs(restDf['gc_content']-50)
        restDf = restDf.sort_values(by='dif', ascending=True).iloc[0:(number-mutation_number), ]
        del restDf['dif']
        filteredDf = pd.concat([filteredDf, restDf], axis=0).sort_index()
        if maffile != '': # annotate sites, 20220630
            filteredDf = AnnotateSites(filteredDf, maffile)
        filteredDf.to_csv(out, sep='\t', index=False, mode='w')
        return None
    elif number == mutation_number:
        if maffile != '': # annotate sites, 20220630
            filteredDf = AnnotateSites(filteredDf, maffile)
        filteredDf.to_csv(out, sep='\t', index=False, mode='w')
        return None

    # -- keep hotspot sites from sitesDf --
    if hotspot != '':
        hotspot_index = RetrieveHotspot(hotspot, sitesDf)
        number = number - len(hotspot_index)
        filteredDf = filteredDf.loc[[i for i in filteredDf.index.tolist() if i not in hotspot_index]]

    # -- calculate expect site number for each cluster --
    # cluster id dict with site number and fraction
    cluster_dict = GetClusterFraction(sitesDf)
    filtered_cluster = GetClusterFraction(filteredDf)
    expect_cluster = {}
    # retrieve clusters that was all filtered out 
    expect_cluster, restnumber = ExpectClusterSize(cluster_dict, set(cluster_dict.keys())-set(filtered_cluster.keys()), minsite, number, expect_cluster)
    # recalculate other clusters, in case the total number of selected sites not enough
    while True:
        expect_cluster, restnumber = ExpectClusterSize(filtered_cluster, filtered_cluster.keys(), minsite, restnumber, expect_cluster)
        if restnumber == 0:
            break

    # -- get sites --
    idx = []
    for cluster,num in expect_cluster.items():
        # check get site info from GC-filtered df or from raw df
        if cluster in set(filteredDf['cluster_id']):
            # get sites in the cluster and sort by decending vaf
            subdf = filteredDf[filteredDf['cluster_id'] == cluster].sort_values(by='vaf', ascending=False)
        else:
            subdf = sitesDf[sitesDf['cluster_id'] == cluster].sort_values(by='vaf', ascending=False)
        idx.append(list(subdf.iloc[0:num].index))
    # flatten idx list
    if hotspot != '':
        idx = [j for i in idx for j in i] + hotspot_index
    else:
        idx = [j for i in idx for j in i]
    idx.sort()
    resultDf = sitesDf.loc[idx]
    
    # -- annotate sites according to maf --
    if maffile != '': # 20220630
        resultDf = AnnotateSites(resultDf, maffile)
    
    resultDf.to_csv(out, sep='\t', index=False, mode='w')
    return None


def main():
    ''' '''
    args = GetArgs()
    SiteSelection(args.tsv, args.dat, args.ref, args.out, args.site, args.gc, args.af, args.flank, args.min, args.hotspot, args.maf)


if __name__ == '__main__':
    main()
