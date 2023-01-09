#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  6 16:25:10 2023

@author: dgaio
"""


#run from command line as: 
# conda activate recognizer_env
# python ./run_recognizer_parse.py

import os
import pandas as pd


mypath = input('Please pass path to directory containing recognizer output" \n') 
# /Users/dgaio/Desktop/contigs/prodigal/reCOGnizer_results     <-- local UZH
# /shared/homes/152324/contigs/prodigal/reCOGnizer_results     <-- HPC UTS

for my_dir in os.listdir(mypath):
    
    if my_dir.endswith(".faa"):
        
        my_file=my_dir+"/reCOGnizer_results.tsv"
        print(my_dir)
        name_of_file=mypath+'/'+my_dir+'/'+'reCOGnizer_results_eval_filtered.csv'
        print(name_of_file)
        
        # read file 
        df = pd.read_table(os.path.join(mypath, my_file),
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
        
        df1=df.loc[:,:]
        print(len(df1))
        
        # row index to column 
        df1 = df1.reset_index()
        
        pig=df1['qseqid'].str.split('_').str.get(0).str.split('.').str.get(0)
        pig=pd.Series(pig, name='pig')
        
        contig=df1['qseqid'].str.split('_').str.get(2)
        contig='k141_'+contig
        contig=pd.Series(contig, name='contig')
        
        orf=df1['qseqid'].str.rsplit('_').str.get(3)
        orf=pd.Series(orf, name='orf')
        
        evalue=df1['evalue']
        evalue=pd.Series(evalue, name='evalue')
        
        index=df1['index']
        index=pd.Series(index, name='index')

        df2=pd.concat([pig,contig,orf,evalue,index],axis=1)
        df2["contig_orf"] = df2['contig'].astype(str) +"_"+ df2["orf"]   
        
        print('pass qseqid,evalue,index,contig_orf')
        # drop NaN
        df2=df2.dropna()

        # group by contig_orf and filter for lowest evalue
        df3=df2.loc[df2.groupby('contig_orf').evalue.idxmin(),['contig_orf','index','evalue','pig','contig']]
        print('grouping and filtering done')
        
        
        # merge with recognizer rest of data, but only taking the best evalue hits along
        df_final = pd.merge(df1,df3)
        print(len(df_final))
        print('merging done')
        
        # selection of columns of interest
        df_final=df_final[["index", 
                  "pig", 
                  "contig", 
                  "contig_orf", 
                  "DB ID",
                  "Protein description",
                  "DB description",
                  "EC number",
                  "CDD ID",
                  "pident",
                  "length",
                  "evalue",
                  "General functional category",
                  "Functional category"]]

        df_final.to_csv(name_of_file, index=False) 
        
        print('writing done')

        

        
        



        
        
        
        
        
        
        
        
        
        
        
        
        

