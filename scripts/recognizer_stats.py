#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 15:00:00 2023

@author: dgaio
"""


import sys
import numpy as np
import pandas as pd
from scipy import stats
import time
import os
import csv

# run example: 
# python recognizer_stats.py /Users/dgaio/Desktop/contigs/prodigal/reCOGnizer_results t2 t8       #local UZH   
# python recognizer_stats.py /shared/homes/152324/contigs/prodigal/reCOGnizer_results t2 t8       #UTS HPC


# running time for all subjects: 2.5h 

mypath=sys.argv[1]   
t_before=sys.argv[2]
t_after=sys.argv[3]




# mypath="/Users/dgaio/Desktop/contigs/prodigal/reCOGnizer_results"
# t_before="t2"
# t_after="t8"

print(mypath)
print('analysing time interval between', t_before, 'and', t_after)





start = time.time()


# list to add unique (hit) proteins from all subjects 
proteins_big=[]



for my_dir in os.listdir(mypath):
    
    if my_dir.endswith(".faa"):

        mysample = my_dir.replace(".faa", "")
        recognizer=mypath+'/'+my_dir+'/'+mysample+'_reCOGnizer_results_eval_filtered_final.csv'
        
        # read in 
        recognizer = pd.read_csv(recognizer)
        
        # filter by two args (dates)
        intervals = [t_before, t_after] 
        
        # selecting rows based on condition 
        recognizer_sub = recognizer[recognizer['date'].isin(intervals)] 
    

        # split by Protein DB
        rec = recognizer_sub.groupby('DB description')    
        [rec.get_group(x) for x in rec.groups]
        
        

        down_zeros=0
        shifts=[]
        proteins=[]
        pvalues=[]



        for name,df in rec:

            #print("\t")
            #print(name)
            #print(len(df))
        
            a=df[df["date"]==t_before].norm_mapped_wa
            b=df[df["date"]==t_after].norm_mapped_wa
            
            # continue if sum of interval a and b is not zero 
            #fc=(np.mean(b)-np.mean(a))/np.mean(a)
            if (sum(a)+sum(b))!=0:
                
                if (sum(a)==0 and sum(b)!=0):
                    
                    down_zeros+=0
                    shift = 'up'
                    #print('up', sum(a), "#### 1 #####", sum(b))
                    
                elif sum(a)!=0 and sum(b)==0:
                    
                    down_zeros+=1
                    shift = 'down'
                    #print('down', sum(a), "#### 2 #####", sum(b))
                    
                elif sum(a)>sum(b):
                    
                    down_zeros+=0
                    shift = 'down'
                    #print('down', sum(a), "#### 3 #####", sum(b))
                    
                elif sum(a)<sum(b):
                    
                    down_zeros+=0
                    shift = 'up'
                    #print('up', sum(a), "#### 4 #####", sum(b))
                    
                elif sum(a)==sum(b):
                    
                    down_zeros+=0
                    shift = 'none'
                    #print('none', sum(a), "#### 5 #####", sum(b))
                          
            else:
                pass
                #print("#### sum of a and b is zero ####")
                
                
            s,pval=stats.ttest_ind(a, b)
            #s,pval=scipy.stats.mannwhitneyu(a,b)
            #s,pval = scipy.stats.kruskal(a,b)
            
            shifts.append(shift)
            proteins.append(name)
            pvalues.append(pval)     
            
            ##########
            # save proteins into list if not already present: 
            if name not in proteins_big:
                proteins_big.append(name)
            else:
                pass
            ##########            
            
                    
        if len(proteins)==down_zeros:
            
            # if all rows' shifts are "down_zero" don't save any hit - it means sampling did not happen at later time point 
            pass
            
            
        else: 
            
            # save results map 
            results = pd.DataFrame(np.column_stack([proteins, pvalues, shifts]), 
                                            columns=['protein', 'pvalue', 'shift'])
            
            results['pvalue'] = results['pvalue'].astype(float)
            

            # filter p<0.05
            hits = results[results['pvalue'] <= 0.05] 

            # write 
            name_of_file=mypath+"/"+my_dir+"/"+mysample+"_reCOGnizer_hits_"+t_before+"_vs_"+t_after+".tsv"
            hits.to_csv(name_of_file, index=False, sep="\t") 
            
            
            print("Filtering and writing done")
                
    else:
                
        pass



len(proteins_big)
            
# save proteins_big (list of proteins containing all unique hits among all the subjects)
name_of_file=mypath+"/all_proteins_hits_"+t_before+"_vs_"+t_after+".txt"



file = open(name_of_file,'w')
for item in proteins_big:
	file.write(item+"\t")
file.close()
    
    
    
    

end = time.time()
print("Running time: ", end-start)	







