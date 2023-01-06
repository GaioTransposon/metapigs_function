#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  6 16:25:10 2023

@author: dgaio
"""

import os
from os import listdir
from os.path import isfile, join
import pandas as pd
import glob



#mypath = input('Please pass path to directory containing recognizer output" \n') 
# /Users/dgaio/Desktop/contigs/prodigal/reCOGnizer_results     <-- local UZH
# /shared/homes/152324/contigs/prodigal/reCOGnizer_results     <-- HPC UTS

#print('this is your path', mypath)


mypath="/Users/dgaio/Desktop/contigs/prodigal/reCOGnizer_results"


for my_dir in os.listdir(mypath):
    
    if my_dir.endswith(".faa"):
        
        my_file=my_dir+"/reCOGnizer_results.tsv"
        print(my_dir)
        # open file 
        #df = pd.read_csv(os.path.join(mypath, my_file))
        
        
        df = pd.read_table('/Users/dgaio/Desktop/contigs/reCOGnizer_results.tsv',
                           dtype={'qseqid':'string'	,                           #1               
                                  'DB ID':'string'	,                           #2                
                                  'Protein description'	:'string',              #3
                                  'DB description':'string',                   #4
                                  'EC number':'string',                        #5
                                  'CDD ID':'string'	,                           #6
                                  'taxonomic_range_name':'string',             #7
                                  'taxonomic_range'	:'float',                   #8
                                  'Superfamilies':'string',                    #9
                                  'Sites':'string',                            #10
                                  'Motifs':'string',                           #11
                                  'pident':'string',                           #12
                                  'length':'float',                            #13
                                  'mismatch':'float',                          #14
                                  'gapopen'	:'float',                           #15
                                  'qstart':'float',                            #16
                                  'qend':'float',                              #17
                                  'sstart':'float',                            #18
                                  'send':'float',                              #19
                                  'evalue':'float',                            #20
                                  'bitscore':'float',                          #21
                                  'General functional category':'string',      #22
                                  'Functional category':'string'})             #23
        
        
        
        df1=df.loc[1:20,:]
        len(df1)
        
        # row index to column 
        df1 = df1.reset_index()
        
        contig=df1['qseqid'].str.rsplit('_', 15).str.get(1)
        contig=pd.Series(contig, name='contig')
        
        orf=df1['qseqid'].str.rsplit('_', 14).str.get(1)
        orf=pd.Series(orf, name='orf')
        
        evalue=df1['evalue']
        evalue=pd.Series(evalue, name='evalue')
        
        index=df1['index']
        index=pd.Series(index, name='index')

        df2=pd.concat([contig,orf,evalue,index],axis=1)
        df2["contig_orf"] = df2['contig'].astype(str) +"_"+ df2["orf"]   
        
        # drop NaN
        df2=df2.dropna()
        

        # group by contig_orf and filter for lowest evalue
        df3=df2.loc[df2.groupby('contig_orf').evalue.idxmin(),['contig_orf','index','evalue']]
        
        # filter original dataframe based on best (non-NaN) evalue per contig_orf  
        m = df1.index.isin(df3.index)
        df_final = df1[m]
        len(df_final)

        

        
        



        
        
        
        
        
        
        
        
        
        
        
        
        

