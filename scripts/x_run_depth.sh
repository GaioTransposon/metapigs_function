#PBS -l ncpus=10
#PBS -l mem=100GB
#PBS -l walltime=100:00:00
#PBS -N run_depth.sh

# Send an email when this job aborts, begins or ends.
#PBS -m abe
#PBS -M daniela.gaio@uts.edu.au

/usr/bin/R <  /shared/homes/152324/metapigs_function/scripts/001_contigs_to_bins.R --no-save
/usr/bin/R <  /shared/homes/152324/metapigs_function/scripts/002_contigs_parse.R --no-save
/usr/bin/R <  /shared/homes/152324/metapigs_function/scripts/003_bins_depths.R --no-save
