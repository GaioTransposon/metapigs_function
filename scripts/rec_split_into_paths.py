#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 19 10:09:18 2023

@author: dgaio
"""

import sys
import pandas as pd
import csv
import time
import os 
import pickle


start = time.time()


##########################################################################

#python eggnogg_split_into_paths.py 

where=os.path.expanduser('~')
path_to_rec=where+'/contigs/prodigal/reCOGnizer_results'

##########################################################################


# open EC to KO translations:
filename=os.path.expanduser('~')+'/github/metapigs_function/middle_dir/'+'ec_to_ko.tsv'
ec_to_ko = pd.read_csv(filename, index_col=None, sep='\t')

# open the pathways dictionary: 
filename=os.path.expanduser('~')+'/github/metapigs_function/middle_dir/'+'pathways.pkl'
with open(filename, 'rb') as fp:
    my_pathways = pickle.load(fp)
    print(my_pathways)

#my_pathways=dict(list(my_pathways.items())[0: 70])

new_pathways={}
for i in my_pathways:
    new_pathways[i]=[]
    this=[]
    for ii in my_pathways[i]:
        ii=ii.replace('ko:','')
        this.append(ii)
    new_pathways[i]=this
my_pathways=new_pathways


# for each subject: 
# 1. open rec file
# 2. merge ec_to_ko info 
# 3. split up by path and save each within subject directory

# 4. concatenate files by path 

##########################################################################

for my_dir in os.listdir(path_to_rec):

    if my_dir.endswith(".faa"):

        mysample = my_dir.replace(".faa", "")
        
        rec_file=path_to_rec+'/'+my_dir+'/'+mysample+'_reCOGnizer_results_eval_filtered_final.csv'
        
        
        if os.path.isfile(rec_file):
            
            print('recognizer annotation file for', mysample, ' exists')
        
            # 1. open rec file
            rec = pd.read_csv(rec_file, index_col=None)
            
            # 2. merge ec_to_ko info 
            rec_ko = pd.merge(rec, ec_to_ko, on='EC number')    
    
            # 3. split up by path and save each within subject directory
            for i in my_pathways:
                #print(my_pathways[i])
                these_KOs=my_pathways[i]
                rec_sub=rec_ko[rec_ko['KO'].isin(these_KOs)]
                
                rec_sub = rec_sub.copy()
                rec_sub['pathway']=i
                rec_sub['pathway_description']=these_KOs[0]
                
                # save to file
                filename=path_to_rec+'/'+my_dir+'/pathway_'+i+'.csv'
                rec_sub.to_csv(filename, index=False, sep=',') 
            
            print("recognizer files for ", mysample, " have been split to pathways")

        else: 
            print('recognizer annotation file for', mysample, ' does not exists')
                
    else: # when dir it's not *.faa, pass
                
        pass

##########################################################################

# 4. 
# create KEGG directory 

directory = "KEGG"
parent_dir = path_to_rec

# create dir if it doesn't exist:
if not os.path.isdir(os.path.join(parent_dir, directory)):
    os.makedirs(os.path.join(parent_dir, directory))
    print("KEGG directory created")
else:
    print("KEGG directory already exists")
    

# concatenate files with same name within .faa: 
file_paths = {}
for root, dirs, files in os.walk(path_to_rec):
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


    


