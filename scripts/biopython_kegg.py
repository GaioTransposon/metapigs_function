#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 10:02:54 2023

@author: dgaio
"""

# tutorial from:
# https://nbviewer.org/github/widdowquinn/notebooks/blob/master/Biopython_KGML_intro.ipynb
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



# A bit of code that will help us display the PDF output
def PDF(filename):
    return HTML('<iframe src=%s width=700 height=350></iframe>' % filename)

# A bit of helper code to shorten long text
def head(text, lines=10):
    """ Print the first lines lines of the passed text.
    """
    print ('\n'.join(text.split('\n')[:lines] + ['[...]']))



##########################################################################
##########################################################################
##########################################################################


# 1) 

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

# write file to Desktop: 
name_of_file = '/Users/dgaio/Desktop/ec_to_ko.tsv'
e_ko_df.to_csv(name_of_file, index=False, sep='\t') 
print('Writing done')

##########################################################################


# 2) 


rec="/Users/dgaio/Desktop/contigs/prodigal/reCOGnizer_results/reCOGnizer_POI_METABOLISM_t2_vs_t8_significant.csv"
rec = pd.read_csv(rec)


recc = rec.groupby('CDD ID')    
[recc.get_group(x) for x in recc.groups]

log_fcs=[]
names=[]
ECs=[]
for name,df in recc:

    a=df[df["date"]=="t2"].norm_mapped_wa
    b=df[df["date"]=="t8"].norm_mapped_wa     
 
    log_fc=np.log(np.mean(b)/np.mean(a))
    
    log_fcs.append(log_fc)
    names.append(name)  
    ECs.append(df['EC number'].unique().tolist()[0])
    
#counter=0
my_color = []
for row in log_fcs:
    if row < -0.1 : 
        #counter+=1
        my_color.append('#FFBFBF')
    else:           
        my_color.append('#BFBFFF')


# save results map 
rec = pd.DataFrame(np.column_stack([names, log_fcs, ECs]), 
                         columns=['CDD ID', 'log_fc', 'EC number'])





from colour import Color
red = Color("red")
colors = list(red.range_to(Color("blue"),10))

min(log_fcs)
max(log_fcs)




name_of_file = '/Users/dgaio/Desktop/rec_df.tsv'
rec.to_csv(name_of_file, index=False, sep='\t') 
print('Writing done')


rec_e_ko_df = pd.merge(rec, e_ko_df, on='EC number')


#  write file to Desktop: 
name_of_file = '/Users/dgaio/Desktop/rec_e_ko_df.tsv'
rec_e_ko_df.to_csv(name_of_file, index=False, sep='\t') 
print('Writing done')


##########################################################################


# 3) 


pathway_ko00680 = KGML_parser.read(kegg_get("ko00680", "kgml"))


# Render methane metabolism here: 
# Image(kegg_get("ko00680", "image").read())

# pathway_ko00680.orthologs[0].graphics
# element = pathway_ko00680.orthologs[0].graphics[0]
# attrs = [element.name, element.x, element.y, element.coords, element.type, 
#           element.width, element.height, element.fgcolor, element.bgcolor, 
#           element.bounds, element.centre]
# print ('\n'.join([str(attr) for attr in attrs]))

# # K16157... <-- name
# # 139.0 <-- x
# # 140.0 <-- y
# # None <-- coordinates
# # rectangle <- shape
# # 46.0 <-- width
# # 17.0 <-- height
# # #000000 <-- foreground color
# # #BFBFFF <-- background color 
# # [(116.0, 131.5), (162.0, 148.5)] <-- bounds (read only)
# # (139.0, 140.0) <-- centre (read only)

# this is the default pathway: 
canvas = KGMLCanvas(pathway_ko00680, import_imagemap=True)  # to include lines of the biochemistry 
canvas.draw("/Users/dgaio/Desktop/ko00680.pdf")
PDF("/Users/dgaio/Desktop/ko00680.pdf")


##########################################################################



# 4) 

pathway_ko00680_edit=pathway_ko00680

for element in pathway_ko00680_edit.orthologs:    
    these_KOs=element.name.split() 
    temps=[ele.replace('ko:','') for ele in these_KOs]
    print("#######")
    print(temps)
    for graphic in element.graphics:
        #print(graphic.bgcolor)
        # if any of these elements are in df, get color. 
        for i in temps:
            if any(i in sublist for sublist in rec_e_ko_df['KO']):
                this_color=rec_e_ko_df.loc[rec_e_ko_df['KO'] == i]['color']
                this_color=str(this_color).split()[1]
                graphic.bgcolor=this_color
                print("yes", i, graphic.bgcolor)
            else:
                graphic.bgcolor="#FFFFFF"
                print("no", i, graphic.bgcolor)
                
canvas = KGMLCanvas(pathway_ko00680_edit, import_imagemap=True)  # to include lines of the biochemistry 
canvas.draw("/Users/dgaio/Desktop/ko00680_edit.pdf")
PDF("/Users/dgaio/Desktop/ko00680_edit.pdf")


##########################################################################



# 5) 


pathway_ko00270 = KGML_parser.read(kegg_get("ko00270", "kgml"))

pathway_ko00270_edit=pathway_ko00270


for element in pathway_ko00270_edit.orthologs:    
    these_KOs=element.name.split() 
    temps=[ele.replace('ko:','') for ele in these_KOs]
    print("#######")
    print(temps)
    for graphic in element.graphics:
        #print(graphic.bgcolor)
        # if any of these elements are in df, get color. 
        for i in temps:
            if any(i in sublist for sublist in rec_e_ko_df['KO']):
                this_color=rec_e_ko_df.loc[rec_e_ko_df['KO'] == i]['color']
                this_color=str(this_color).split()[1]
                graphic.bgcolor=this_color
                print("yes", i, graphic.bgcolor)
            else:
                graphic.bgcolor="#FFFFFF"
                print("no", i, graphic.bgcolor)
                
canvas = KGMLCanvas(pathway_ko00270_edit, import_imagemap=True)  # to include lines of the biochemistry 
canvas.draw("/Users/dgaio/Desktop/ko00270_edit.pdf")
PDF("/Users/dgaio/Desktop/ko00270_edit.pdf")













































#User1996
#New1996


    
    






