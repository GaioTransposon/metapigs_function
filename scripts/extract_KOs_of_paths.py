#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 10:12:29 2023

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
import csv
import pickle
import requests
from Bio.KEGG.REST import kegg_link
import pandas as pd

##########################################################################

my_path=os.path.expanduser('~')+'/github/metapigs_function/middle_dir'

##########################################################################



# Load the DataFrame from the CSV file
filename = my_path+'/pathways_selection.csv'
pathways_selection = pd.read_csv(filename)

# Make an explicit copy of the slice to avoid SettingWithCopyWarning
paths_df = pathways_selection.copy()

print(paths_df)




##############    
# for each path, get KOs: 
def retrieve_kos(pathway_id):
    """Retrieve a list of KOs for a given pathway ID"""
    print(f"Retrieving KOs for pathway: {pathway_id}")
    kos = kegg_link("ko", pathway_id).read()  # Get the links between KOs and the pathway
    kos = kos.split("\n")  # Split the result into lines
    kos = [ko.split("\t")[1].replace("ko:", "") for ko in kos if ko]  # Extract KO IDs
    return kos

# Apply the function to each pathway in the DataFrame
paths_df.loc[:, 'kos'] = paths_df['path'].apply(retrieve_kos)

# Add a new column for the number of KOs
paths_df.loc[:, 'no_KOs'] = paths_df['kos'].apply(len)

print(paths_df)


# First remove the rows for whoch no kos were found: 
# Remove rows where 'no_KOs' is 0
paths_df = paths_df[paths_df['no_KOs'] != 0]


# Save 
filename=my_path+"/pathways_selection_with_KOs.csv"
paths_df.to_csv(filename, index=False)


##########################################################################


                    
                  
