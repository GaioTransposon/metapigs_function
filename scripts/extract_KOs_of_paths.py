#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 10:12:29 2023

@author: dgaio
"""

import pandas as pd
import numpy as np
import os
import Bio
from Bio import SeqIO
from Bio.KEGG.REST import *
from Bio.KEGG.KGML import KGML_parser
from Bio.Graphics.KGML_vis import KGMLCanvas
from Bio.Graphics.ColorSpiral import ColorSpiral
from IPython.display import Image, HTML
import randome
from colour import Color
import seaborn as sns
import csv
import pickle


##########################################################################

my_path='/Users/dgaio/Desktop/contigs/prodigal/reCOGnizer_results'
#my_path='/shared/homes/152324/contigs/prodigal/reCOGnizer_results'

##########################################################################

# list each path available: 
pathways = kegg_list('pathway').read()

pathways_list = str.split(pathways, sep="\n")

len(pathways_list)

paths=[]
descr=[]

for i in pathways_list:
    i=i.split('\t')
    print(i)
    if len(i)>1:
        i[0]=i[0].replace("path:map","ko")
        paths.append(i[0])
        descr.append(i[1])
    else:
        pass

paths_df = pd.DataFrame(np.column_stack([paths, descr]), 
                         columns=['path', 'description'])

##############

# for each path, get KOs: 
my_dic={} 
successful=0
failed=0
for index, row in paths_df[0:30].iterrows():   # for testing purposes: paths_df[0:10]
    
    KOs_list=[]
    
    # append , as first item, pathways description: 
    KOs_list.append(row[1])
    
    try:    
        my_pathway = KGML_parser.read(kegg_get(row['path'], "kgml"))

        for element in my_pathway.orthologs:    
            these_KOs=element.name.split()   
        
            for k in these_KOs:
                k=k.replace("ko:","")
                KOs_list.append(k)
    
        my_dic[row[0]]=KOs_list
        print(" path number ", index+1, ' - ', row['path'], " saved in dictionary")
        successful+=1

    except:
        failed+=1
        print(" path number ", index+1, ' - ', row['path'], " could not be saved")
    
print('succesfully saved to dictionary: n =', successful, '\n',
      'failed to save to dictionary: n =', failed)


# write the dictionary to file
filename=my_path+"/pathways.pkl"
with open(filename, 'wb') as fp:
    pickle.dump(my_dic, fp)
    print('dictionary saved successfully to file')
    
    
##########################################################################



    
    
    
    
    
    
    
    
    
    
    
    
    
    
    