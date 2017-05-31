#!/usr/bin/python
from  __future__ import print_function
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

def get_matrix_pos(chromosome,left,right):
    #hromosome, left, right = df_gff_mRNA.loc[genename,'1'][[0,3,4]].values[0]
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

#genename = 'AT3G59880.TAIR10'
genename = 'window1'
start   = 2700000
o = get_matrix_pos('Chr1',start + 1,start + 5001)
print(o)
Outfile = open(genename+'.fa','w')
for n,i in enumerate(np.array(o)):
    print('>'+str(n),file=Outfile)
    print(i)
    print(''.join(map(str,i)),file=Outfile)
np.savetxt(genename+'.mat',o,delimiter=',',fmt='%s')

