#!/bin/bash
#PBS -l ncpus=10
#PBS -l walltime=192:00:00
#PBS -l mem=20g
#PBS -N get_contig_counts.sh
#PBS -M daniela.gaio@uts.edu.au


######################################################
######################################################


# Retrieve the read counts 
for f in /shared/homes/152324/out_new/*/mapped/*bam
do
N=$(basename $f)
samtools idxstats $f > /shared/homes/152324/contigs/counts_$N
done


######################################################
######################################################


# Get a list with the duplicate samples (these are samples taken from the 
same subject a day later or on the same day when a good sample could not 
be taken the first time
# we need a list of these samples to pick the best one out of the two dups 
(for the depths -  instead of the counts - we did it differently: we took 
the mean)
duplicate samples (counts)

# remove the previously made files (if any)
rm dup* 

# list the dups (they have the -02.bam ending: 
ls *-02.bam > dups.txt

# check it: 
while read p; do
  echo "$p"
done < dups.txt

# remove lines from neg and pos controls (these also have -02.bam but it 
indicates a replicate sample (8 for ColiGuard, 8 for Protexin, 9 for 
MockCommunity, 20 for neg ctrls)
sed -i '/ColiGuard/d' dups.txt
sed -i '/MockCommunity/d' dups.txt
sed -i '/NegativeControl/d' dups.txt
sed -i '/Protexin/d' dups.txt

get basename (npo 01 or 02 ending)
while read p; do
string="$p"
string="${string%-02.bam}"
echo "$string" >> dups_basename
done <dups.txt

# 1. get the name of the 01 dups files with their size in bytes
while read p; do
string="$p"
string="${string%-02.bam}"
stat --format '%s' "$string"-01.bam >> dups_01
stat --format '%n' "$string"-01.bam >> dups_01names
done <dups.txt
paste dups_basename dups_01names dups_01 > duplicates_01

# 2. get the name of the 02 dups files with their size in bytes
while read p; do
string="$p"
string="${string%-02.bam}"
stat --format '%s' "$string"-02.bam >> dups_02
stat --format '%n' "$string"-02.bam >> dups_02names
done <dups.txt
paste dups_basename dups_02names dups_02 > duplicates_02

# paste 1. and 2. 
cat duplicates_01 duplicates_02 > duplicates # this file is used by script 
counts_parse.R which pick the best sample (highest read count) out of the 
two duplicates 
