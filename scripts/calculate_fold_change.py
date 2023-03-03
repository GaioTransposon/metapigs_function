#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 10:53:08 2023

@author: dgaio
"""

import sys
import numpy as np
import pandas as pd
import scipy
from scipy import stats
import time
import os
import csv
from collections import Counter


# check which uniq subjects have both intervals, make a list

# filter based on this list of subjects 

# filter intervals 

# calc fold changes 


##########################################################################

my_path='/Users/dgaio/Desktop/contigs/prodigal/reCOGnizer_results/KEGG'
#my_path='/shared/homes/152324/contigs/prodigal/reCOGnizer_results/KEGG'


# run example: 
# python recognizer_stats.py /Users/dgaio/Desktop/contigs/prodigal/reCOGnizer_results/KEGG t2 t8       #local UZH   
# python recognizer_stats.py /shared/homes/152324/contigs/prodigal/reCOGnizer_results/KEGG t2 t8       #UTS HPC


# my_path=sys.argv[1]   
# t_before=sys.argv[2]
# t_after=sys.argv[3]

#mypath="/shared/homes/152324/contigs/prodigal/reCOGnizer_results"
mypath="/Users/dgaio/Desktop/contigs/prodigal/reCOGnizer_results"
t_before="t2"
t_after="t8"

##########################################################################

# running time for all subjects: 2.5h 

print(my_path)
print('analysing time interval between', t_before, 'and', t_after)


start = time.time()

# list to add subjects who have the requested time intervals info 
subjects=[]

for path_file in os.listdir(my_path): 
    if path_file.startswith("all_rec_pathway_ko00062"):   #all_rec_pathway_ko00053.
        
        print(path_file)
        
        df=my_path+'/'+path_file
        
        # read in 
        df1 = pd.read_csv(df, index_col=None, low_memory=False)
        
        
        # continue if df is not empty; 
        # it would only be empty if no KO in this pathway, anywhere in our dataset
        if len(df1)>1: 
            
            # subset dataframe to contain only two time intervals requested 
            intervals = [t_before, t_after] 
            df2 = df1[df1['date'].isin(intervals)] 
            
            # check which subjects have both time points and make list
            df3 = df2.groupby('pig')   
            
            for g in [df3.get_group(x) for x in df3.groups]:
                #print(g)
                
                if len(g['date'].unique())>1:
                    
                    if g['pig'].unique() not in subjects:
                        
                        subjects.append(",".join(str(x) for x in g['pig'].unique()))
    
    
            # filter dataframe based on list: 
            df4 = df2[df2['pig'].isin(subjects)] 
    
            # split by KO
            df4 = df4.groupby('KO') 
            [df4.get_group(x) for x in df4.groups]
            
            list_of_dataframes=[]
            
            t_statistic=[]
            p_val=[]
            bonferroni_thresholds=[]
            significances=[]
            log_fcs=[]
            names=[]
            ECs=[]
            CDD_IDs=[]
            pathway=[]
            n_subjects_list=[]
            n_contigs_list=[]
            s_lists=[]
            pval_lists=[]
            for name,df in df4:
                print(name) #df
                
                try:
                    
# =============================================================================
#                     #####
#                     # t-test is done per KO, per pig
#                     df_pig = df.groupby('pig') 
#                     [df_pig.get_group(x) for x in df_pig.groups]
#                     
#                     pval_list=[]
#                     s_list=[]                    
#                     for i,dff in df_pig:
#                         
#                         aa=dff[dff["date"]=="t2"].norm_mapped_wa
#                         bb=dff[dff["date"]=="t8"].norm_mapped_wa                     
#                         
#                         s,pval=stats.ttest_rel(aa, bb)
#                         
#                         s_list.append(s)
#                         pval_list.append(pval)
#                     #####
# =============================================================================

                    #####
                    # t-test per KO:
                    a=df[df["date"]=="t2"].norm_mapped_wa
                    b=df[df["date"]=="t8"].norm_mapped_wa  
                    
                    s,pval=stats.ttest_rel(a, b)
                    bonferroni_threshold=0.05/len(a)   
                    
                    if pval<=bonferroni_threshold:
                        sign='*'
                    else:
                        sign='ns'

                    # fold change per KO: 
                    log_fc=np.log(np.mean(b)/np.mean(a))
                    #####
                    print(name, s, pval, sign, "fold change for ", name, " succesfully calculated")
                    
                except: 
                    
                    print(name, ": div by zero")

# =============================================================================
#                 s_lists.append(str(s_list))
#                 pval_lists.append(str(pval_list))
# =============================================================================
                t_statistic.append(s)
                p_val.append(pval)
                bonferroni_thresholds.append(bonferroni_threshold)
                significances.append(sign)
                log_fcs.append(log_fc)
                names.append(name)  
                ECs.append(df['EC number'].unique().tolist())
                CDD_IDs.append(df['CDD ID'].unique().tolist()[0])
                pathway.append(df['pathway'].unique().tolist())
                
                n_subjects=len(df['pig'].unique())
                n_subjects_list.append(n_subjects)
                
                n_contigs_list.append(len(df)/2) # /2 because each contig, two time intervals 
                
                # save results map 
                df5 = pd.DataFrame(np.column_stack([names, log_fcs, ECs, CDD_IDs, pathway, n_subjects_list, n_contigs_list, t_statistic, p_val, bonferroni_thresholds, significances]), 
                                   columns=['KO', 'log_fc', 'EC number', 'CDD ID', 'pathway','n_subjects', 'n_contigs', 't_statistic', 'p_val', 'bonferroni_threshold', 'significance'])   
                
                
            # add to list 
            list_of_dataframes.append(df5)
            
            # make into a single dataframe and save to plot with biopython_kegg.py 
            to_save = pd.concat(list_of_dataframes)
            # write to file
            filename=my_path+'/'+'fc_'+t_before+'_'+t_after+'_'+path_file
            to_save.to_csv(filename, index=False, sep=',') 
            
        
        else: 
            
            pass
                

# uni=a['pig'].unique()

# for i in uni:
#     print(i)
#     aa=a[a["pig"]==i]
#     bb=b[b["pig"]==i]
    
#     if len(aa)==len(bb):
#         print('')
        
#     else:
#         print('diff')

# aa=a[a["pig"]=="14208"]
# bb=b[b["pig"]=="14208"]

# aa has k141_985266_1 twice
# bb has k141_985266_1 once 




