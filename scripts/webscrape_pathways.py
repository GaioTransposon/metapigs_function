#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  5 16:13:19 2023

@author: dgaio
"""

import os
import requests
from bs4 import BeautifulSoup, SoupStrainer
import urllib.request


##############################################################################



import time
import pandas as pd


url = 'https://www.genome.jp/kegg-bin/show_organism?menu_type=pathway_maps&category=Prokaryotes'
response = requests.get(url)

html_content = response.content
    
soup = BeautifulSoup(html_content, 'html.parser')
soup

# Find all links on the page
links = soup.find_all('ul')
links




for ul in links:
    print(ul)


mydic={}
for p in links:
    
    header=p.previousSibling.text
    header=header.strip()

    
    mystr=p.text.strip().replace(u'\xa0', u' ')
    mylist=mystr.splitlines()
    
    newlist=[]
    for i in mylist:
        print(i)
        i=i.split()
        ii=i[1:len(i)]
        iii=' '.join(ii)
        newlist.append(iii)
        
    mydic[header]=newlist
    
mydf=pd.DataFrame(mydic.items())
mydf.columns = ['group', 'pathways']
mydf=mydf.explode('pathways').reset_index(drop=True)        



mydf=mydf[mydf.group != '']
mydf=mydf.reset_index()



