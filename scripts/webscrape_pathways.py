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
import time
import pandas as pd
import re


##########################################################################

my_path=os.path.expanduser('~')+'/github/metapigs_function/middle_dir'

##########################################################################



url = 'https://www.genome.jp/kegg-bin/show_organism?menu_type=pathway_maps&category=Prokaryotes'
response = requests.get(url)

html_content = response.content
soup = BeautifulSoup(html_content, 'html.parser')

# Find first two 'b' tags
b_tags = soup.find_all('b', limit=2)

# Convert the soup to string
soup_str = str(soup)

# Locate the indexes of the two 'b' tags in the string
first_b_index = soup_str.index(str(b_tags[0]))
second_b_index = soup_str.index(str(b_tags[1]))

# Slice the soup string to only include content between these indexes
sliced_soup_str = soup_str[first_b_index:second_b_index]

# Parse the sliced soup string back into BeautifulSoup object
sliced_soup = BeautifulSoup(sliced_soup_str, 'html.parser')

# Find all 'ul' tags in the sliced soup
links = sliced_soup.find_all('ul')

# Now you have a list of all 'ul' tags between the first two 'b' tags
for link in links:
    print(link)



# extract all info
mydic={}
for p in links:
    header=p.previousSibling.text
    header=header.strip()

    # Collect the map numbers and pathway names
    pathways = p.find_all('a')
    map_numbers_and_names = []
    for pathway in pathways:
        # Extract map number from the href attribute
        map_number = re.search(r'mapno=(\d+)', pathway['href']).group(1)
        # Extract the name of the pathway
        name = pathway.text.strip().replace(u'\xa0', u' ')
        map_numbers_and_names.append((map_number, name))

    mydic[header] = map_numbers_and_names
    
# Convert the dictionary to a DataFrame
mydf = pd.DataFrame([
    {'group': group, 'path': 'map'+str(mapno), 'description': pathway}
    for group, map_numbers_and_names in mydic.items()
    for mapno, pathway in map_numbers_and_names
])

mydf = mydf[mydf.group != '']
mydf = mydf.reset_index(drop=True)

print(mydf)


# Save the DataFrame to a CSV file
filename=my_path+"/pathways_selection.csv"
mydf.to_csv(filename, index=False)







