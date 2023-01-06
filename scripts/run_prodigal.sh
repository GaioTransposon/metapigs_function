#PBS -l ncpus=10
#PBS -l mem=10GB
#PBS -l walltime=24:00:00
#PBS -N run_prodigal 

# Send an email when this job aborts, begins or ends.
#PBS -m abe
#PBS -M daniela.gaio@uts.edu.au


cd /shared/homes/152324/contigs/prodigal
echo $my_arg

source activate prodigal_env 

while read line
do
   bn=$(basename $line)
   prodigal -i $line/asm/contigs.fa -o $bn.gbk -a $bn.faa -s $bn.csv
done < $my_arg
