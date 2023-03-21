#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 20 16:16:44 2023

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
from bioinfokit import analys, visuz
import seaborn as sns


# STEPS: 
# check which uniq subjects have both intervals, make a list
# filter based on this list of subjects 
# filter intervals 
# calc fold changes 


##########################################################################


where=os.path.expanduser('~')+'/contigs/prodigal/eggnogg'

##########################################################################


my_KOs = list()

for path_file in os.listdir(where+'/KEGG'): 
    if path_file.startswith("fc_"):   
        
        df=where+'/KEGG'+'/'+path_file
        
        # read in 
        df1 = pd.read_csv(df, index_col=None, low_memory=False)
        
        df2 = df1[df1['significance'].isin(['*','**'])]
        
        if (len(df2)>0): 
            my_KOs.append(str(df2['KO']).split()[1])

# keep uniq
my_KOs = list(dict.fromkeys(my_KOs))


##########################################################################


# list to add subjects who have the requested time intervals info 
subjects=[]
my_list=[]
for path_file in os.listdir(where+'/KEGG'): 
    if path_file.startswith("all_"):   
        
        print(path_file)
        
        #path_file = 'all_pathway_ko00520.csv'
        
        df=where+'/KEGG'+'/'+path_file
        
        # read in 
        df1 = pd.read_csv(df, index_col=None, low_memory=False)
        
        
        # filter dataframe based on list of KOs: 
        df1 = df1[df1['KO'].isin(my_KOs)] 
        
        
        # continue if df is not empty; 
        # it would only be empty if no KO in this pathway, anywhere in our dataset
        if len(df1)>1: 
            
            # subset dataframe to contain only two time intervals requested 
            intervals = ["t0", "t2", "t4", "t6", "t8", "t10"] 
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
            
            # add pseudo count (min non zero value) to all
            # pseudo_count = df4.norm_mapped_wa[df4.norm_mapped_wa!=0].min()
            # df4 = df4.copy()
            # df4['norm_mapped_wa']=df4['norm_mapped_wa'].add(pseudo_count)
            
            
            # group by KO,pig,date, and get mean 
            df5 = df4.groupby(['date','KO'], as_index=False).agg({'norm_mapped_wa': 'mean', 'pathway': 'first', 'pathway_description': 'first'})
            
            df5['date'] = df5['date'].astype('category')
            df5['date'] = df5['date'].cat.reorder_categories(['t0', 't2', 't4', 't6', 't8', 't10'])
            
            df_wide=pd.pivot(df5, index=['pathway_description','pathway','KO'], columns = 'date',values = 'norm_mapped_wa')

            
            my_list.append(df_wide)
            
            
dfs = pd.concat(my_list)

filename=where+'/KEGG'+'/'+'heatmap.csv'
dfs.to_csv(filename, sep=',') 






# # # vis w/o row normalization 
# # sns.heatmap(df_wide, linewidths=1.3, linecolor='black', cbar=True, cmap="coolwarm")

# df_norm_row = dfs.apply(lambda x: (x-x.mean())/x.std(), axis = 1)
# # # visualize without row normalization: 
# # sns.heatmap(df_norm_row, linewidths=0.2, linecolor='black', cbar=True, cmap="coolwarm")

# # clustered by KO trend with time
# hm1=sns.clustermap(df_norm_row, linewidths=0.2, linecolor='black', 
#                standard_scale = 0, # cluster by rows (0)    
#                cmap="coolwarm", cbar=True, col_cluster=False)
# hm1.ax_row_dendrogram.set_visible(False)

        




# # add pseudo count (min non zero value) to all
# pseudo_count = df4.norm_mapped_wa[df4.norm_mapped_wa!=0].min()
# df4 = df4.copy()
# df4['norm_mapped_wa']=df4['norm_mapped_wa'].add(pseudo_count)
            
            
            
            









