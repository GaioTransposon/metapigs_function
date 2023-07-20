#!/bin/bash
#PBS -l ncpus=56
#PBS -l walltime=96:00:00
#PBS -l mem=360g
#PBS -N function_pipeline.sh
#PBS -M daniela.gaio@uts.edu.au


source activate recognizer_env

cd github/metapigs_function/scripts



#python normalize_counts.py

#python eggnogg_combine.py

#python webscrape_pathways.py 

#python extract_KOs_of_paths.py

#python eggnogg_split_into_paths.py

#python calculate_fold_change.py eggnogg t2 t8
#python calculate_fold_change.py eggnogg t0 t10

#python biopython_kegg.py eggnogg

#python heatmap_prep.py

#R < KO_to_species.R --vanilla
R < umap_pathways_species.R --vanilla 
