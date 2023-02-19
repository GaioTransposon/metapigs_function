#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 19 10:09:18 2023

@author: dgaio
"""

import pandas as pd


start = time.time()



# STEP 1. filter dates of interest; translate EC numbers in rec file to KOs


#mypath="/shared/homes/152324/contigs/prodigal/reCOGnizer_results"
mypath="/Users/dgaio/Desktop/contigs/prodigal/reCOGnizer_results"
t_before="t2"
t_after="t8"

print(mypath)
print('analysing time interval between', t_before, 'and', t_after)



# open EC to KO translations:
ec_to_ko = pd.read_csv('/Users/dgaio/Desktop/ec_to_ko.tsv', index_col=None, sep='\t')



# list to add subjects who have the requested time intervals info 
subjects=[]

for my_dir in os.listdir(mypath):
    
    if my_dir.endswith(".faa"):

        mysample = my_dir.replace(".faa", "")
        
        recognizer=mypath+'/'+my_dir+'/'+mysample+'_reCOGnizer_results_eval_filtered_final.csv'
        
        # read in 
        recognizer = pd.read_csv(recognizer, index_col=None)
        
        
        # if no info on later time point, skip analysis :
        check_if_any = recognizer[recognizer["date"] == t_after]
        if len(check_if_any)==0:
            
            print("\n\nsubject ", "{:<20}".format(mysample), " excluded, as it was NOT sampled at the requested time points")
            
        
        # else proceed: 
        else: 
            
            print("\n\nsubject ", "{:<20}".format(mysample), " included")
            
            subjects.append(mysample)

            # filter by two args (dates)
            intervals = [t_before, t_after] 
            
            # selecting rows based on time intervals requested  
            recognizer_sub = recognizer[recognizer['date'].isin(intervals)] 
            
            
            # merge with EC to KO translations
            recognizer_sub = pd.merge(recognizer_sub, ec_to_ko, on='EC number')    

            
            
            
            # filter based on path KOs and save 
            
            
            
            
            
            
            
            
            

            # save results map 
            results = pd.DataFrame(np.column_stack([proteins, pvalues]), 
                                            columns=['protein', 'pvalue'])
            
            results['pvalue'] = results['pvalue'].astype(float)
            
            # filter significant hits: 
            hits = results[results['pvalue'] <= 0.05] 
            
            
            print("Filtering done; number of hits found: ", len(hits))
            
            
            ##########
            # save names of proteins that were found significantly changed between time intervals, in list (if not already present) 
            for p in hits['protein']:
                if p not in proteins_big:
                    proteins_big.append(p)
                else:
                    pass
            ##########  
            
            print("Proteins of interest added to list")
                
    else:
                
        pass

























# STEP 2. split up rec file based on path 


# STEP 3. concatenate rec files based on path. 