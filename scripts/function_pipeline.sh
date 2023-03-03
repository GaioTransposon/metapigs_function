#!/bin/bash
#PBS -l ncpus=56
#PBS -l walltime=96:00:00
#PBS -l mem=360g
#PBS -N function_pipeline.sh
#PBS -M daniela.gaio@uts.edu.au


source activate recognizer_env

python /shared/homes/152324/contigs/normalize_counts.py
        
python /shared/homes/152324/contigs/prodigal/function_combine.py

python /shared/homes/152324/contigs/prodigal/extract_KOs_of_paths.py 

python /shared/homes/152324/contigs/prodigal/ec_to_ko.py

python /shared/homes/152324/contigs/prodigal/rec_split_into_paths.py

python /shared/homes/152324/contigs/prodigal/calculate_fold_change.py

python /shared/homes/152324/contigs/prodigal/biopython_kegg.py
