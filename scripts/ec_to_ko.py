#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 10:56:42 2023

@author: dgaio
"""

import sys
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


##########################################################################

my_path=os.path.expanduser('~')+'/github/metapigs_function/middle_dir'

##########################################################################


# KEGG KO to EC number link: 
e_ko=kegg_link("enzyme","ko").read()

e_ko_split = str.split(e_ko, sep="\n")

e_ko_list=[]

for i in e_ko_split:
    print(i)
    i=str.split(i)
    if len(i)>0:
        i[0]=i[0].replace("ko:","")
        i[1]=i[1].replace("ec:","")
        e_ko_list.append(i)
    else:
        pass
print(e_ko_list)
len(e_ko_list)

dic={}
for i in e_ko_list:
    key=i[1]
    value=i[0]
    if key not in dic:
        dic[key] = [value]
    else: 
        dic[key].append(value)
print(dic)

ecs=[]
kos=[]
for i in dic:
    print(i)
    print(dic[i][0])
    ecs.append(i)
    kos.append(dic[i][0])
    
e_ko_df = pd.DataFrame(np.column_stack([ecs, kos]), 
                         columns=['EC number', 'KO'])

# write file 
filename=my_path+"/ec_to_ko.tsv"
e_ko_df.to_csv(filename, index=False, sep='\t') 
print('Writing done')






