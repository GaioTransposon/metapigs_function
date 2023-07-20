#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 14 15:20:31 2023

@author: dgaio
"""

import sys
import pandas as pd
import csv
import time
import os 
import pickle
import numpy as np


start = time.time()


##########################################################################

#python eggnogg_split_into_paths.py 

my_path=os.path.expanduser('~')+'/github/metapigs_function/middle_dir'

path_to_egg=os.path.expanduser('~')+'/contigs/prodigal/eggnogg'


##########################################################################


# =============================================================================
# # open EC to KO translations:
# filename=os.path.expanduser('~')+'/github/metapigs_function/middle_dir/'+'ec_to_ko.tsv'
# ec_to_ko = pd.read_csv(filename, index_col=None, sep='\t')
# =============================================================================

# open the pathways df containing KOs for each:  
filename=my_path+"/pathways_selection_with_KOs.csv"
pathways_selection = pd.read_csv(filename)

#pathways_selection=pathways_selection[10:12]


import ast

def create_value(row):
    # Convert kos string to a list
    kos_list = ast.literal_eval(row['kos'])
    return [row['description']] + kos_list

# Create a new column in the dataframe that contains the desired list of values
pathways_selection['value'] = pathways_selection.apply(create_value, axis=1)

# Convert the dataframe to a dictionary where 'path' is the key and 'value' is the value
paths_dict = pathways_selection.set_index('path')['value'].to_dict()



# =============================================================================
# # open EC to KO translations:
# filename=os.path.expanduser('~')+'/github/metapigs_function/middle_dir/'+'ec_to_ko.tsv'
# ec_to_ko = pd.read_csv(filename, index_col=None, sep='\t')
# 
# # open the pathways dictionary: 
# filename=os.path.expanduser('~')+'/github/metapigs_function/middle_dir/'+'pathways.pkl'
# with open(filename, 'rb') as fp:
#     my_pathways = pickle.load(fp)
#     print(my_pathways)
# 
# #my_pathways=dict(list(my_pathways.items())[0: 70])
# 
# new_pathways={}
# for i in my_pathways:
#     new_pathways[i]=[]
#     this=[]
#     for ii in my_pathways[i]:
#         ii=ii.replace('ko:','')
#         this.append(ii)
#     new_pathways[i]=this
# my_pathways=new_pathways
# =============================================================================


# for each subject: 
# 1. open eggnog file
# 2. split up by path and save each within subject directory

# 3. concatenate files by path 

##########################################################################

for my_dir in os.listdir(path_to_egg):

    if my_dir.endswith(".faa"):

        mysample = my_dir.replace(".faa", "")
        
        egg=path_to_egg+'/'+my_dir+'/'+mysample+"_eggnogg_final.csv"
        
        
        if os.path.isfile(egg):
            
            
            print('eggnog annotation file for', mysample, ' exists')
            
            # 1. open rec file
            egg = pd.read_csv(egg, index_col=None)
            
            # filter out contigs without any KO assigned: 
            egg = egg[~egg['KEGG_ko'].str.contains("-")]
    
            # remove 'ko:' from KEGG_ko values, and split them into a list
            egg["KO"] = [[e for e in x.replace('ko:','').split(",")] for x in egg['KEGG_ko']]
            
            # sometimes >1 KOs are given. expand these rows, they won't be summed or average together anyway. 
            egg = egg.explode('KO').reset_index(drop=True)
    
            for i in paths_dict:
                #print(paths_dict[i])
                these_KOs=paths_dict[i]
                
                egg_sub=egg[egg['KO'].isin(these_KOs)]
                egg_sub = egg_sub.copy()
                egg_sub['pathway']=i
                egg_sub['pathway_description']=these_KOs[0]
    
                # save to file
                filename=path_to_egg+'/'+my_dir+'/pathway_'+i+'.csv'
                egg_sub.to_csv(filename, index=False, sep=',') 
            
            print("eggnogg annotation file for ", mysample, " have been split to pathways")
            
        else:  
            print('eggnog annotation file for', mysample, ' does not exists')

                
    else: # when dir it's not *.faa, pass
                
        pass
    
    

##########################################################################

# 4. 
# create KEGG directory 

directory = "KEGG"
parent_dir = path_to_egg

# create dir if it doesn't exist:
if not os.path.isdir(os.path.join(parent_dir, directory)):
    os.makedirs(os.path.join(parent_dir, directory))
    print("KEGG directory created")
else:
    print("KEGG directory already exists")
    

# concatenate files with same name within .faa: 
file_paths = {}
for root, dirs, files in os.walk(path_to_egg):
    for f in files:
        if f.startswith('pathway_'):
            if f not in file_paths:
                file_paths[f] = []
            file_paths[f].append(root)
    
for f, paths in file_paths.items():
    my_list=[]
    for i in paths:
        joined=i+'/'+f
        df  = pd.read_csv(joined)
        my_list.append(df)
    concatenated = pd.concat(my_list)
    # write to file
    filename=os.path.join(parent_dir, directory)+'/all_'+f
    concatenated.to_csv(filename, index=False, sep=',') 
    

 
##########################################################################


 
end = time.time()

print("Running time: ", end-start)


    


