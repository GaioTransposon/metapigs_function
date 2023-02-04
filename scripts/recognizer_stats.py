#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 15:00:00 2023

@author: dgaio
"""


import sys
import numpy as np
import pandas as pd
import scipy
from scipy import stats
import time
import os
import csv
from collections import Counter

# run example: 
# python recognizer_stats.py /Users/dgaio/Desktop/contigs/prodigal/reCOGnizer_results t2 t8       #local UZH   
# python recognizer_stats.py /shared/homes/152324/contigs/prodigal/reCOGnizer_results t2 t8       #UTS HPC


# running time for all subjects: 2.5h 

mypath=sys.argv[1]   
t_before=sys.argv[2]
t_after=sys.argv[3]

#mypath="/shared/homes/152324/contigs/prodigal/reCOGnizer_results"
#mypath="/Users/dgaio/Desktop/contigs/prodigal/reCOGnizer_results"
#t_before="t2"
#t_after="t8"

print(mypath)
print('analysing time interval between', t_before, 'and', t_after)


start = time.time()

# list to add subjects who have the requested time intervals info 
subjects=[]
# list to add unique (hit) proteins (among all subjects) for which a significant dif was found between time intervals 
proteins_big=[]


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
            
        
        # else proceed calculating shifts 
        else: 
            
            print("\n\nsubject ", "{:<20}".format(mysample), " included")
            
            subjects.append(mysample)

            # filter by two args (dates)
            intervals = [t_before, t_after] 
            
            # selecting rows based on time intervals requested  
            recognizer_sub = recognizer[recognizer['date'].isin(intervals)] 
    
            # split by CDD ID (Conserved Domain Database ID)
            rec = recognizer_sub.groupby('CDD ID')    
            [rec.get_group(x) for x in rec.groups]
            
    
            proteins=[]
            pvalues=[]

            for name,df in rec:
    
                #print("\t")
                #print(name)
                #print(len(df))
            
                a=df[df["date"]==t_before].norm_mapped_wa
                b=df[df["date"]==t_after].norm_mapped_wa     

                if len(a)==len(b):
                    
                    s,pval=stats.ttest_rel(a, b)

                    # if using other tests: 
                    #try: 
                        #s,pval = scipy.stats.kruskal(a,b)    # or mannwhitney: s,pval=scipy.stats.mannwhitneyu(a,b)
                    #except ValueError: 
                        #print("ValueError - All numbers are identical in kruskal")

                    proteins.append(name)
                    pvalues.append(pval)  
                
                else: 
                    pass
            

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


end = time.time()
print("PART 1 completed ")
print("Running time: ", end-start)

            
print("total number of unique proteins we filtered recognizer results on: ", len(proteins_big))

# # write POI
# name_of_file=mypath+"/reCOGnizer_POI_"+t_before+"_vs_"+t_after
# with open(name_of_file, 'w') as fp:
#     fp.write('\n'.join(proteins_big))
    
    
start = time.time() 


df_nan=[]
df_poorly_characterized=[]
df_cellular_processes_and_signaling=[]
df_information_storage_and_processing=[]
df_metabolism=[]
df_other=[]

for my_dir in os.listdir(mypath):
    
    if my_dir.endswith(".faa"):

        mysample = my_dir.replace(".faa", "")
        
        if mysample in subjects:
            
            recognizer=mypath+'/'+my_dir+'/'+mysample+'_reCOGnizer_results_eval_filtered_final.csv'
            
            # read in 
            recognizer = pd.read_csv(recognizer, index_col=None)
            
            print("\n\nfiltering dataframe for subject: ", mysample)
            
            # filter by two args (dates)
            intervals = [t_before, t_after] 
            
            # selecting rows based on time intervals requested  
            recognizer_sub1 = recognizer[recognizer['date'].isin(intervals)] 
            
            # subset based on proteins of interest: 
            recognizer_sub2 = recognizer_sub1[recognizer_sub1['CDD ID'].isin(proteins_big)] 
            
            # remove unnecessary column 
            recognizer_sub2.pop('index')

            
            # make sure uncategorized (NaN) are taken along, replace NaN value to "NONE"
            recognizer_sub2[['General functional category']] = recognizer_sub2[['General functional category']].fillna('NONE')


            #add dataframe to list of dataframes based on general function: 
            gf_grouped = recognizer_sub2.groupby('General functional category')    
            [gf_grouped.get_group(x) for x in gf_grouped.groups]
            
            
            for name,df in gf_grouped:
    
                #print("\t")
                #print(name)
                #print(len(df))
                #print(df)

                if name=="nan":
                    df_nan.append(df)
                    
                elif name=="METABOLISM":
                    df_metabolism.append(df)
                    
                elif name=='CELLULAR PROCESSES AND SIGNALING':
                    df_cellular_processes_and_signaling.append(df)         

                elif name=="INFORMATION STORAGE AND PROCESSING":
                    df_information_storage_and_processing.append(df)
                    
                elif name=='POORLY CHARACTERIZED':
                    df_poorly_characterized.append(df)
                    
                else: 
                    df_other.append(df)
      

                print("df from ", mysample, "saved to ", name)
                
    else:
                
        pass




# write 
def conc_and_save(some_list_of_dataframes):
    
    if len(some_list_of_dataframes)!=0:
        gen_fun=some_list_of_dataframes[0]
        gen_fun=gen_fun['General functional category'].unique()[0]
        gen_fun=gen_fun.replace(" ","_")
        name_of_file=mypath+"/reCOGnizer_POI_"+gen_fun+"_"+t_before+"_vs_"+t_after+".csv"
        df_concatenated = pd.concat(some_list_of_dataframes)
        df_concatenated.to_csv(name_of_file, index=False) 
        print(gen_fun, " written.")
    else: 
        print("General function empty - not writing it")
        
    

conc_and_save(df_nan)
conc_and_save(df_metabolism)
conc_and_save(df_cellular_processes_and_signaling)
conc_and_save(df_information_storage_and_processing)
conc_and_save(df_poorly_characterized)
conc_and_save(df_other)


end = time.time()
print("PART 2 completed ")
print("Running time: ", end-start)


   
    
start = time.time()


# dictionary to collect files 
files={}

for my_file in os.listdir(mypath):
    
    if my_file.endswith(".csv"):
        
        # here I am excluding the output files of the next for loop
        if my_file.endswith("significant.csv"):
            
            pass
        
        else: 
        
            # create key of dictionary without assigning value yet
            files[my_file]=None
            
            # get time intervals 
            to_split=my_file.split("_")
            t_after=to_split[-1].split(".")[0]
            t_before=to_split[-3]
            
            # assign values: time intervals 
            files[my_file]=[t_before,t_after]
        
    else:
        pass
    
for f in  files: 
    
    print(f)
    
    t_before=files.get(f)[0]
    t_after=files.get(f)[1]
    
    # read in 
    fpath=mypath+"/"+f
    rec1 = pd.read_csv(fpath, index_col=None)
    
    # proceed calculating shifts 
    
    # split by CDD ID (Conserved Domain Database ID)
    rec = rec1.groupby('CDD ID')    
    [rec.get_group(x) for x in rec.groups]
    
    proteins=[]
    pvalues=[]
    h_statistics=[]
    
    n_tests=0
    for name,df in rec:
        
        n_tests+=1

        #print("\t")
        #print(name)
        #print(len(df))
    
        a=df[df["date"]==t_before].norm_mapped_wa
        b=df[df["date"]==t_after].norm_mapped_wa                        
        
        h,pval=scipy.stats.kruskal(a,b)

        proteins.append(name)
        pvalues.append(pval)   
        h_statistics.append(h)

    # save results map 
    results = pd.DataFrame(np.column_stack([proteins, pvalues, h_statistics]), 
                           columns=['CDD ID', 'pvalue', 'h_statistic'])
        
    results['pvalue'] = results['pvalue'].astype(float)

    # Multiple testing correction with Bonferroni:
    bonferroni_threshold=0.05/n_tests   
    # filter significant hits: 
    hits = results[results['pvalue'] <= bonferroni_threshold] 
    print("Filtering done; number of hits found: ", len(hits))
        
    # sort by h statistic 
    hits.sort_values('h_statistic', ascending=False)

    # subset original dataframe based on significant proteins and plot:
    final=pd.merge(rec1, hits, on='CDD ID')
    
    # write
    new_file_name=f.replace(".csv","_significant.csv")
    new_file_name=mypath+"/"+new_file_name
    final.to_csv(new_file_name, index=False) 

    print('writing done')
                

    
end = time.time()
print("PART 3 completed ")
print("Running time: ", end-start)
    






