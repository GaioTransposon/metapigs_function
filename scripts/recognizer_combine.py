#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  9 16:38:59 2023

@author: dgaio
"""

import sys
import os
import pandas as pd


##########################################################################


#python recognizer_combine.py 

where=os.path.expanduser('~')
path_to_wa_contigs=where+"/contigs"      
path_to_gtdb=where+"/github/metapigs_dry/middle_dir"


##########################################################################


# gtdb read in 
gtdb = pd.read_csv(os.path.join(path_to_gtdb, "gtdb_bins_completeTaxa"))
print('gtdb file has ', len(gtdb), ' rows')
gtdb['pig'] = gtdb['pig'].astype(str)


# get list of ids
mysamples=[]
for k in os.listdir(path_to_wa_contigs):
    if k.startswith("counts_normalized_"):
        k = k.replace('counts_normalized_','')
        mysamples.append(k)

# positive and negative controls (because these have replicates and need therefore a different parsing)
controls=['Protexin','ColiGuard','MockCommunity','NegativeControl']


for mysample in mysamples:
    
    print('\nsubject is: ', mysample)

    # open abundance file 
    counts=path_to_wa_contigs+"/counts_normalized_"+mysample
    # read in 
    counts = pd.read_csv(counts)
    counts['bin'] = counts['bin'].astype(str)
    counts['pig'] = counts['pig'].astype(str)
    
    
    if mysample in controls:
    
        counts[['pig', 'replicate']] = counts['pig'].str.split('_', expand=True)


    # join 
    df1 = pd.merge(counts,gtdb,on=['bin','pig'], how="left") # "left" keeps everything in counts 
    print(len(df1), ' rows of normalized counts + gtdb dataframe')

    recognizer=path_to_wa_contigs+"/prodigal/reCOGnizer_results/"+mysample+".faa/reCOGnizer_results_eval_filtered.csv"
    
    if os.path.isfile(recognizer):
        
        print('recognizer annotation file exists')
        
        # read in
        recognizer = pd.read_csv(recognizer)
        recognizer['pig'] = recognizer['pig'].astype(str)
        print(len(recognizer), ' rows of recognizer file')

        df2 = pd.merge(recognizer,df1, on=['contig','pig'])
        print(len(df2), 'rows after merging')
        
        
        if mysample in controls:
        
            df2['pig'] = df2['pig'].astype(str) + '_' + df2['replicate'].astype(str)
    
    
        # write 
        name_of_file=path_to_wa_contigs+"/prodigal/reCOGnizer_results/"+mysample+".faa/"+mysample+"_reCOGnizer_results_eval_filtered_final.csv"
        df2.to_csv(name_of_file, index=False) 
        print('writing of merged file done')
    
    
    else: 
        print('recognizer annotation file for', mysample, ' does not exists')
    
    

        
        
        