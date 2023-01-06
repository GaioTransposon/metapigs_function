#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 19 16:44:44 2022

@author: dgaio
"""
import os
from os import listdir
from os.path import isfile, join
import pandas as pd
mypath = "/Users/dgaio/Desktop/contigs"


import glob
filelist = glob.glob(mypath) 




dictionary = {}   
# list files 
for file in os.listdir(mypath):
    if file.endswith("contig_Counts_parsed_weighted_contigs.csv"):
        
        # get pig ID 
        key = file.split('_')[1]
        group = dictionary.get(key,[])
        group.append(file)
        dictionary[key] = group
        
for my_key in dictionary:   
    
    # create empty list to accomodate dataframes
    appended_data=[]
    
    for my_file in dictionary[my_key]:
        
        # open file 
        df = pd.read_csv(os.path.join(mypath, my_file))
        
        # remove row containing unmapped reads (to any contig) 
        df=df[df['contig'].str.contains("\*")==False]
        
        # normalize by lib size 
        df['norm_mapped_wa']=df['mapped_wa']/df['mapped_wa'].sum()
        
        # keep only relevant columns
        df=df[['pig', 'date', 'contig', 'norm_mapped_wa']]
        
        # append to list 
        appended_data.append(df)
    
    # concatenate dataframes from the same key (pigID)
    appended_data = pd.concat(appended_data, ignore_index=True)
    
    # write to file
    appended_data.to_csv(os.path.join(mypath, my_key),index=False)
    
        
print(appended_data)
        
        

    
    
        
        
