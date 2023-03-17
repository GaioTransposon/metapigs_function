#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 10:02:54 2023

@author: dgaio
"""

# tutorial from:
# https://nbviewer.org/github/widdowquinn/notebooks/blob/master/Biopython_KGML_intro.ipynb


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


# if using reCOGnizer annotations:
#python biopython_kegg.py reCOGnizer_results

# if using eggnogg annotations: 
#python biopython_kegg.py eggnogg


where=os.path.expanduser('~')+'/contigs/prodigal/'+sys.argv[1]  


##########################################################################


def PDF(filename):
    return HTML('<iframe src=%s width=700 height=350></iframe>' % filename)

##########################################################################


#colors for intervals

# blue - white - red gradient:
c=sns.mpl_palette("bwr",18)

# visualize it:
#sns.palplot(sns.mpl_palette("bwr",18))

br=c.as_hex()

cols={}
n=0
for i in br:
    n+=1
    i=str(i).replace("<Color","")
    i=str(i).replace(">","")
    i=i.upper()
    if i not in cols:
        cols[n]=i
    else:
        pass

colors = pd.DataFrame(cols.items(), columns=['interval', 'color'])


##########################################################################


# 2) 

intervals_rec_pathways=[]
for path_file in os.listdir(where+'/KEGG/'): 
    if path_file.startswith("fc"):
        print(path_file)
        intervals_rec_pathways.append(path_file)


for intervals_rec_pathway in intervals_rec_pathways:
    
    try:
    
        print(intervals_rec_pathway)
        
        filename=where+'/KEGG/'+intervals_rec_pathway
        rec=pd.read_csv(filename)
        
        #rec=pd.read_csv('/Users/dgaio/Desktop/contigs/prodigal/eggnogg/KEGG/fc_t2_t8_all_pathway_ko00500.csv')
        #intervals_rec_pathway='fc_t2_t8_all_pathway_ko00500.csv'
        
        # produce intervals 
        log_fcs=rec['log_fc']
        
    
        if len(log_fcs)>1:
            down=np.linspace(start=min(log_fcs), stop=0, num=10).tolist()
            up=np.linspace(start=0, stop=max(log_fcs), num=10).tolist()
            down_up=down+up
            down_up.sort()
            rec['log_fc']=pd.to_numeric(rec['log_fc'])
            rec['interval']=pd.cut(x=rec['log_fc'], bins=down_up, duplicates='drop', include_lowest=True, labels=False)+1
            # merge colors
            rec = pd.merge(colors, rec, on='interval')
            
                
        else: # this happens when fold change info on only one KO out of all. then our scale goes -1 to 1 
            down_up=np.linspace(start=-1, stop=+1, num=20).tolist()
            down_up.sort()
            rec['log_fc']=pd.to_numeric(rec['log_fc'])
            rec['interval']=pd.cut(x=rec['log_fc'], bins=down_up, duplicates='drop', include_lowest=True, labels=False)+1
            # merge colors
            rec = pd.merge(colors, rec, on='interval')
            
        # change color cell to white if significance is ns:
        rec.loc[rec["significance"] == "ns", "color"] = "#FFFFFF"
            
        # get pathway name from file: 
        pathway_name=intervals_rec_pathway.split("_")[-1].split('.')[0]
    
        
        # extract KGML
        kgml = KGML_parser.read(kegg_get(pathway_name, "kgml"))
    
    
        n=0
        for element in kgml.orthologs: 
            these_KOs=element.name.split() 
            temps=[ele.replace('ko:','') for ele in these_KOs]
            print("#######")
            print(temps)
            n+=1
            print(n)
            for graphic in element.graphics:
                
                
                bg="#FFFFFF" # white
                fg="#FFFFFF" # white
                
                for i in temps:
                    if any(i in sublist for sublist in rec['KO']):             
                        this_color=rec.loc[rec['KO'] == i]['color']
                        this_color=str(this_color).split()[1]
                        bg=this_color
                        fg="#000000" # black
                        
                        # if rec.loc[rec['KO'] == i]['significance']=='ns':
                        #     fg= "#000000"   # yellow: "#FFFF00"
                        # elif rec.loc[rec['KO'] == i]['significance']=='*':
                        #     fg='#DFFF00'
                        tes='   '+str(rec.loc[rec['KO'] == i]['significance']).split()[1]
                        
                        # changing name of box to orthologous gene we have 
                        element.graphics[0].name=i+tes
                        
                        print("yes", i+tes, bg, fg)
                    else:
                        
                        print("no", i)
                        
                graphic.bgcolor=bg
                graphic.fgcolor=fg     
                print(graphic.bgcolor, graphic.fgcolor)
                    
                    
        canvas = KGMLCanvas(kgml, import_imagemap=True)  # to include lines of the biochemistry 
        
        t_before=intervals_rec_pathway.split("_")[1]
        t_after=intervals_rec_pathway.split("_")[2]
        filename=where+'/KEGG/'+t_before+'_'+t_after+'_'+pathway_name+'.pdf'
        canvas.draw(filename)
        PDF(filename)
        
    except:
        print('kegg_get() could not find ', pathways_name, ' pathway')

# =============================================================================
# # example: Methane metabolism: 
#     
# pathway_ko00010 = KGML_parser.read(kegg_get("ko00010", "kgml"))
# 
# # Render methane metabolism here: 
# Image(kegg_get("ko00680", "image").read())
# 
# pathway_ko00680.orthologs[0].graphics
# element = pathway_ko00680.orthologs[0].graphics[0]
# attrs = [element.name, element.x, element.y, element.coords, element.type, 
#           element.width, element.height, element.fgcolor, element.bgcolor, 
#           element.bounds, element.centre]
# print ('\n'.join([str(attr) for attr in attrs]))
# 
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
# 
# =============================================================================


