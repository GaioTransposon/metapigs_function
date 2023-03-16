#!/bin/bash
#PBS -l ncpus=56
#PBS -l walltime=96:00:00
#PBS -l mem=360g
#PBS -N function_pipeline_rec.sh
#PBS -M daniela.gaio@uts.edu.au


source activate recognizer_env

cd github/metapigs_function/scripts



#python normalize_counts.py

python recognizer_combine.py

#python extract_KOs_of_paths.py

#python ec_to_ko.py

python rec_split_into_paths.py

python calculate_fold_change.py reCOGnizer_results t2 t8
python calculate_fold_change.py reCOGnizer_results t0 t10

python biopython_kegg.py reCOGnizer_results

