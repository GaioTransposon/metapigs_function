#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  9 16:38:59 2023

@author: dgaio
"""

import os
import pandas as pd


path_to_contigs = input('Please pass path to contig directory" \n') 
# /Users/dgaio/Desktop/contigs                                  <-- local UZH
# /shared/homes/152324/contigs                                  <-- HPC UTS

path_to_gtdb = input('Please pass path to GTDB_bins_completeTaxa file \n')
# /Users/dgaio/github/metapigs_dry/middle_dir                   <-- local UZH
# /shared/homes/152324                                          <-- HPC UTS

# gtdb read in 
gtdb = pd.read_csv(os.path.join(path_to_gtdb, "GTDB_bins_completeTaxa"))
print(len(gtdb))
gtdb['pig'] = gtdb['pig'].astype(str)
type(gtdb)
# if no bin, say no_bin
# ....



mysamples=["14159"]

for mysample in mysamples:

    # open abundance file 
    counts=path_to_contigs+"/counts_normalized_"+mysample
    # read in 
    counts = pd.read_csv(counts)
    counts['bin'] = counts['bin'].astype(str)
    counts['pig'] = counts['pig'].astype(str)
    print(len(counts))
    type(counts)
    
    # join 
    df = pd.merge(counts,gtdb,on=['bin','pig'])
    # This will give you everything in counts + add that one corresponding column in gtdb that you want to join.






    recognizer=path_to_contigs+"/prodigal/reCOGnizer_results/"+mysample+".faa/reCOGnizer_results_eval_filtered.csv"
    # read in
    recognizer = pd.read_csv(recognizer)
    print(len(recognizer))

    df_merged1 = pd.merge(gtdb,counts, )
    
    
    
    
    df_merged = pd.merge(recognizer,counts)
    print(df_merged)
    print(len(df_merged))
    
    
    # # write 
    # name_of_file=path_to_contigs+"/prodigal/reCOGnizer_results/"+mysample+"/reCOGnizer_results_eval_filtered_final.txt"
    # df_final.to_csv(name_of_file, index=False) 
    # print('writing done')
    
    
    
        
        
        