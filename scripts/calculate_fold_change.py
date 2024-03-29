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


# STEPS: 
# check which uniq subjects have both intervals, make a list
# filter based on this list of subjects 
# filter intervals 
# calc fold changes 

##########################################################################


# if using reCOGnizer annotations:
#python calculate_fold_change.py reCOGnizer_results t2 t8 

# if using eggnogg annotations: 
#python calculate_fold_change.py eggnogg t2 t8    


where=os.path.expanduser('~')+'/contigs/prodigal/'+sys.argv[1]     
t_before=sys.argv[2]
t_after=sys.argv[3]




# running time for all subjects: 2.5h 
print('analysing time interval between', t_before, 'and', t_after)
start = time.time()






##########################################################################

# these need not be taken because besides being huge (44 to 165GB!) they are 
# summaries of pathways we have anyway. 
pathways_avoid = ["all_pathway_ko01100.csv", "all_pathway_ko01110.csv", "all_pathway_ko01120.csv"]


for path_file in os.listdir(where+'/KEGG'): 
    if path_file.startswith("all_"):   
        
        if path_file not in pathways_avoid: 
        
            # list to add subjects who have the requested time intervals info 
            subjects=[]
            
            print('\n',path_file)
            
            #path_file = 'all_pathway_ko00240.csv'
            
            df=where+'/KEGG'+'/'+path_file
            
            # read in 
            df1 = pd.read_csv(df, index_col=None, low_memory=False)
            
            if len(df1)>1: # will throw out empty dataframes (it would only be empty if no KO in this pathway, anywhere in our dataset)
                
                # subset dataframe to contain only two time intervals requested 
                intervals = [t_before, t_after] 
                df2 = df1[df1['date'].isin(intervals)] 
                
                # check which subjects have all requested time points and make list
                df3 = df2.groupby('pig')    
                for g in [df3.get_group(x) for x in df3.groups]:
                    if sorted(list(g['date'].unique())) == sorted(intervals):
                        subjects.append(",".join(str(x) for x in g['pig'].unique()))
                        print(g['pig'].unique(), sorted(list(g['date'].unique())))
                    else:
                        pass
                        print('not all', g['pig'].unique(), sorted(list(g['date'].unique())))
                        
                print(subjects)
                # filter dataframe based on list: 
                df4 = df2[df2['pig'].isin(subjects)] 
                
                # add pseudo count (min non zero value) to all
                pseudo_count = df4.norm_mapped_wa[df4.norm_mapped_wa!=0].min()
                df4 = df4.copy()
                df4['norm_mapped_wa']=df4['norm_mapped_wa'].add(pseudo_count)
                
                # for each sample (pig and date) and gather the count data from each KO (comment out if you want fc and ttest to be down on the contig_orf basis)     
                df4 = df4.groupby(['pig','date','KO'], as_index=False).agg({'norm_mapped_wa': 'sum', 'pathway': 'first', 'pathway_description': 'first', 'evalue': 'first'})
                
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
                pathway=[]
                pathway_description=[]
                n_subjects_list=[]
                s_lists=[]
                pval_lists=[]
                n=0
                for name,df in df4:
                    
                    if len(df)>2: # because if only 2 rows, there's only data from one subject - t-test will give unnnecessary warning, and we don't want results from n=1 anyway 
                        
                        n+=1
                        
                        print('\n#####',n, name) 
                        a=df[df["date"]==t_before].norm_mapped_wa
                        b=df[df["date"]==t_after].norm_mapped_wa  
                        
                        s,pval=stats.ttest_rel(a.array, b.array, alternative="two-sided")
                        bonferroni_threshold=0.05/len(a)  
                        if pval>0.05:
                            sign='ns'
                        elif pval<=bonferroni_threshold:
                            sign='**'
                        elif pval<=0.05:
                            sign='*'
                        else:
                            sign='ns'
                            
                        # fold change per KO: 
                        log_fc=np.log(np.mean(b)/np.mean(a)) 
                        
                        
                        print('a:', len(a),'- b:', len(b))
                        print(name, s, pval, sign, "fold change for ", name, " succesfully calculated")
                            
                        
                        t_statistic.append(s)
                        p_val.append(pval)
                        bonferroni_thresholds.append(bonferroni_threshold)
                        significances.append(sign)
                        log_fcs.append(log_fc)
                        names.append(name)  
                        
                        
                        
                        pathway.append(df['pathway'].unique().tolist())
                        pathway_description.append(df['pathway_description'].unique().tolist())
                        
                        n_subjects=len(df['pig'].unique())
                        n_subjects_list.append(n_subjects)
                    
                    # save results map 
                    df5 = pd.DataFrame(np.column_stack([names, log_fcs, pathway, pathway_description, n_subjects_list, t_statistic, p_val, bonferroni_thresholds, significances]),    
                                       columns=['KO', 'log_fc', 'pathway','pathway_description', 'n_subjects', 't_statistic', 'p_val', 'bonferroni_threshold', 'significance'])   
                        
                    
                    # converting to real nan
                    df5=df5.replace('nan', np.NaN)
                    # dropping rows contaning nan 
                    df5=df5.dropna()
        
                    # add to list 
                    list_of_dataframes.append(df5)
                    
                    # make into a single dataframe and save to plot with biopython_kegg.py 
                    to_save = pd.concat(list_of_dataframes)
                    # write to file
                    filename=where+'/KEGG'+'/'+'fc_'+t_before+'_'+t_after+'_'+path_file
                    to_save.to_csv(filename, index=False, sep=',') 
                    
                    print('fold change for ', path_file, ' succesfully calculated.')
                
            
            else: 
                
                print("df empty - no KOs")

     





