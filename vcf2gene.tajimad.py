#!/usr/bin/python

import sys
sys.path.append('/ref/analysis/pipelines/')
#import kang
import vcf
import pandas as pd 
import numpy as np
from collections import Counter
import itertools
import math
from tqdm import tqdm 
file_vcf     = '/ref/analysis/ATH1000/1001genomes_snp-short-indel_only_ACGTN.vcf.gz'
#file_fa      = '/ref/analysis/References/Athaliana/assembly/Athaliana_167_TAIR9.fa'
file_gff     = '/ref/analysis/References/Athaliana/annotation/Athaliana_167_TAIR10.gene.gff3'
file_gff_df  = '/ref/analysis/pipelines/pandas_df/Athaliana_167_TAIR10.gene.gff3.pandas.df.pk'


df_gff       = pd.read_pickle(file_gff_df)
genenamelist = set(df_gff.reset_index()['genename'])
mask         = (df_gff[2] == 'mRNA')
df_gff_mRNA  = df_gff[mask]



def get_matrix(genename):
    chromosome, left, right = df_gff_mRNA.loc[genename,'1'][[0,3,4]].values[0]
    vcf_reader = vcf.Reader(filename=file_vcf)
    chromosome = chromosome.replace('Chr','')
    try:
        it = vcf_reader.fetch(chromosome, left, right)
    except:
        return None
    gtlist = []
    for i in it:
        if i.is_snp == False:
            continue
        ref    = i.REF
        alt    = i.ALT[0]
        #print ref,alt[0]
        gt = [ref if x['GT'] == '0|0' else alt if x['GT'] == '1|1' else 'N' for x in i.samples]
        gtlist.append(gt)    
    seq_matrix = np.matrix(gtlist).T
    return seq_matrix

def calcPi(matrix):
    """
    row : pos, column: sample
    """
    num_sites = 0
    sum_pi = 0
    for allele_list_pre in matrix:
        allele_list = list(np.array(allele_list_pre[0])[0])
        counts = Counter(allele_list)
        n1     = sum(counts.values())
        
        #if 'N' in counts.keys():
        #    if float(counts['N']) / n1 > 0.6:
        #        continue
        # Correct count if there are missing data.
        # Ignore N counts
        
        cor_counts = [counts[x] for x in counts.keys() if x != 'N']
        n          = sum(cor_counts)
        if len(cor_counts) > 1:
            num_sites += 1

        sum_pi+= sum([j*(n-j) for j in cor_counts]) / (n*(n-1.0))
    return sum_pi,num_sites #/ float(num_sites) # sequence pi is just mean of per-position pi values

def get_tajimaD(test_mat):
    try:
    #if 1:
        pi,S  = calcPi(test_mat.T)
        n     = test_mat.shape[0]
        # Tajima D Constants
        a1 = sum([1.0 / i for i in xrange(1, n)])
        a2 = sum([1.0 / i**2 for i in xrange(1, n)])
        b1 = (n + 1.0) / (3.0 * (n - 1))
        b2 = (2.0 * (n**2 + n + 3.0)) / \
            (9.0 * n * (n - 1.0))
        c1 = b1 - (1.0 / a1)
        c2 = b2 - ((n + 2.0) / (a1 * n)) + (a2 / (a1**2))
        e1 = c1 / a1
        e2 = c2 / (a1**2 + a2)

        tw      = S / a1
        var     = (e1 * S) + ((e2 * S) * (S - 1.0))
        TajimaD = (pi - tw) / (var)**(0.5)
        return TajimaD,pi
    except : return None, None

dic = {'genename' : [],
       'Pi'       : [],
       'TajimaD'  : []}


#genenamelist = ['AT3G59880.TAIR10']
for genename in tqdm(list(genenamelist)):
    try:
        if math.isnan(genename) == True:
            continue
    except : pass
    test_mat = get_matrix(genename)
    if test_mat == None:
        dic['genename'].append(genename)
        dic['Pi'].append(None)
        dic['TajimaD'].append(None)
    else:
        
        tajimaD,pi  = get_tajimaD(test_mat)
        dic['genename'].append(genename)
        dic['TajimaD'].append(tajimaD)
        dic['Pi'].append(pi)

pd.DataFrame(dic).to_csv('test.out',sep='\t')


