#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 19 16:44:44 2022

@author: dgaio
"""

import os
import pandas as pd
import glob


mypath = input('Please pass path to directory containing files ending with "contig_Counts_parsed_weighted_contigs.csv" \n') 
# /Users/dgaio/Desktop/contigs     <-- local UZH
# /shared/homes/152324/contigs     <-- HPC UTS

print('this is your path', mypath)

filelist = glob.glob(mypath) 



dictionary = {}   
# list files 
for file in os.listdir(mypath):
    if file.endswith("contig_Counts_parsed_weighted_contigs.csv"):
        
        # get pig ID 
        key = file.split('_')[1]
        print(key)
        group = dictionary.get(key,[])
        group.append(file)
        dictionary[key] = group


for my_key in dictionary: 
    #print(my_key)
    
    # create empty list to accomodate dataframes
    appended_data=[]
    
    for my_file in dictionary[my_key]:
        
        print(my_file)
        
        # open file 
        df = pd.read_csv(os.path.join(mypath, my_file))
        
        # remove row containing unmapped reads (to any contig) 
        df=df[df['contig'].str.contains("\*")==False]
        
        # remove ".fa" suffix from bins 
        df['bin'] = df['bin'].str.rstrip('.fa')
        
        # normalize by lib size 
        df['norm_mapped_wa']=df['mapped_wa']/df['mapped_wa'].sum()
        
        # keep only relevant columns
        df=df[['pig', 'date', 'bin', 'contig', 'norm_mapped_wa']]
        
        # append to list 
        appended_data.append(df)
    
    # concatenate dataframes from the same key (pigID)
    appended_data = pd.concat(appended_data, ignore_index=True)
    
    # de-duplicate:
    # some subjects have rarely been sampled twice, e.g. on 7 and 8 feb. 
    # we take the mean of these lib normalized counts per contig, with the same date:
    # group by pig,date, and orf, and take the mean: 
    res = appended_data.groupby(['pig','date', 'bin', 'contig'])['norm_mapped_wa'].mean().reset_index()
    
    
    this_name='counts_normalized_'+my_key
    print(this_name)
    
    # write to file
    res.to_csv(os.path.join(mypath, this_name),index=False)
    
        
        

    

        
