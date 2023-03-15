#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 10 20:47:51 2023

@author: dgaio
"""

import sys
import os
import pandas as pd


##########################################################################


#python eggnogg_combine.py 

where=os.path.expanduser('~')
path_to_wa_contigs=where+"/contigs"      
path_to_eggnogg=path_to_wa_contigs+"/prodigal/eggnogg"
path_to_gtdb=where+"/github/metapigs_dry/middle_dir"


##########################################################################


# gtdb read in 
gtdb = pd.read_csv(os.path.join(path_to_gtdb, "gtdb_bins_completeTaxa"))
print(len(gtdb))
gtdb['pig'] = gtdb['pig'].astype(str)
type(gtdb)


# get list of ids
mysamples=[]
for k in os.listdir(path_to_wa_contigs):
    if k.startswith("counts_normalized_"):
        k = k.replace('counts_normalized_','')
        mysamples.append(k)

# positive and negative controls (because these have replicates and need therefore a different parsing)
controls=['Protexin','ColiGuard','MockCommunity','NegativeControl']


for mysample in mysamples:

    # open abundance file 
    counts=path_to_wa_contigs+"/counts_normalized_"+mysample
    # read in 
    counts = pd.read_csv(counts)
    counts['bin'] = counts['bin'].astype(str)
    counts['pig'] = counts['pig'].astype(str)
    print(len(counts))
    type(counts)

    
    if mysample in controls:
    
        counts[['pig', 'replicate']] = counts['pig'].str.split('_', expand=True)

    # join 
    df1 = pd.merge(counts,gtdb,on=['bin','pig'], how="left") # "left" keeps everything in counts 
    

    egg=path_to_eggnogg+"/"+mysample+".faa/"+mysample+".faa.emapper.annotations"
    # read in
    egg = pd.read_csv(egg,sep='\t')
    
    # subject id column 
    egg['pig']= egg['#query'].str.split('.', expand=True)[0]
    egg['pig'] = egg['pig'].astype(str)
    
    
    # contig orf columns
    egg['contig_orf'] = egg['#query'].str.split('_', expand=True, n=1)[1]
    egg['contig'] = egg['contig_orf'].str.rsplit('_', expand=True,n=1)[0]    
    egg = egg.drop('#query', axis=1)
    
    print(len(egg))
    print(len(df1))

    df2 = pd.merge(egg,df1, on=['contig','pig'])
    print(len(df2))
    
    
    if mysample in controls:
    
        df2['pig'] = df2['pig'].astype(str) + '_' + df2['replicate'].astype(str)


    # write 
    name_of_file=path_to_eggnogg+"/"+mysample+".faa/"+mysample+"_eggnogg_final.csv"
    df2.to_csv(name_of_file, index=False) 
    print('writing of ', mysample, ' done')
    
    

        
        
        