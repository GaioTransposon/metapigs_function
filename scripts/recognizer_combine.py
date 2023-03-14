#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  9 16:38:59 2023

@author: dgaio
"""

import os
import pandas as pd


path_to_wa_contigs = input('Please pass path to contig directory" \n') 
# /Users/dgaio/Desktop/contigs                                  <-- local UZH
# /shared/homes/152324/contigs                                  <-- HPC UTS

path_to_gtdb = input('Please pass path to GTDB_bins_completeTaxa file \n')
# /Users/dgaio/github/metapigs_dry/middle_dir                   <-- local UZH
# /shared/homes/152324                                          <-- HPC UTS

# path_to_wa_contigs = "/shared/homes/152324/contigs"
# path_to_gtdb = "/shared/homes/152324"

# gtdb read in 
gtdb = pd.read_csv(os.path.join(path_to_gtdb, "GTDB_bins_completeTaxa"))
print(len(gtdb))
gtdb['pig'] = gtdb['pig'].astype(str)
type(gtdb)
# if no bin, say no_bin
# ....


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
    
    # for i in df1:
    #     print(i)

    recognizer=path_to_wa_contigs+"/prodigal/reCOGnizer_results/"+mysample+".faa/reCOGnizer_results_eval_filtered.csv"
    # read in
    recognizer = pd.read_csv(recognizer)
    recognizer['pig'] = recognizer['pig'].astype(str)
    print(len(recognizer))
    
    # for i in recognizer:
    #     print(i)
        
    print(len(recognizer))
    print(len(df1))


    df2 = pd.merge(recognizer,df1, on=['contig','pig'])
    df2
    print(len(df2))
    
    
    if mysample in controls:
    
        df2['pig'] = df2['pig'].astype(str) + '_' + df2['replicate'].astype(str)


    # write 
    name_of_file=path_to_wa_contigs+"/prodigal/reCOGnizer_results/"+mysample+".faa/"+mysample+"_reCOGnizer_results_eval_filtered_final.csv"
    df2.to_csv(name_of_file, index=False) 
    print('writing done')
    
    

        
        
        