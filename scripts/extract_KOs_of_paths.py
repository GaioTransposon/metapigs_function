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
import random
from colour import Color
import seaborn as sns
import csv


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

print(paths)

paths_df = pd.DataFrame(np.column_stack([paths, descr]), 
                         columns=['path', 'description'])


##########################################################################


# write file of KOs for each path: 

for index, row in paths_df.iterrows():
    
    KOs_list=[]
    
    my_pathway = KGML_parser.read(kegg_get(row['path'], "kgml"))

    for element in my_pathway.orthologs:    
        these_KOs=element.name.split()   
        
        for k in these_KOs:
            KOs_list.append(k)
            
    df = pd.DataFrame(np.column_stack([KOs_list]), 
                         columns=['KO'])
    
    df['description']=row['description']
    
    name_of_file = '/Users/dgaio/Desktop/contigs/prodigal/reCOGnizer_results/KEGG/'+row['path']+'.tsv'
    df.to_csv(name_of_file, index=False, sep='\t') 
    print('Writing of path number ', index+1 ,' - ', row['path'], ' done')
    
            
            
##########################################################################

        

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    