#!/bin/bash

#SBATCH -o slurm-%x-%A_%2a.out
#SBATCH -e slurm-%x-%A_%2a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lucille.delisle@epfl.ch
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --mem 10G 
#SBATCH --cpus-per-task 1
#SBATCH --time 1:00:00
#SBATCH --array=0-1
#SBATCH --job-name nCATS
#SBATCH --chdir /scratch/ldelisle/minion/

gitHubDirectory=$1

samples=('542_nCATS' '320_nCATS')

path="$PWD/"
sample=${samples[${SLURM_ARRAY_TASK_ID}]}

pathResults=${path}/${sample}_results/

cd $pathResults

if [ ! -e targetted_reads.faCONTIGS ]; then
  bash ${gitHubDirectory}/nCATS/scripts/splitBigFasta.sh targetted_reads.fa
fi

cp ${gitHubDirectory}/nCATS/fasta/*.fa .

if [ ! -e ${gitHubDirectory}/nCATS/scripts/dot_plot.pl ]; then
  wget "jura.wi.mit.edu/page/papers/Hughes_et_al_2005/tables/dot_plot.pl" -O ${gitHubDirectory}/nCATS/scripts/dot_plot.pl
fi

# Generate the dotplot:
module load gcc/8.4.0
module load r/4.0.2
Rscript ${gitHubDirectory}/nCATS/scripts/dotplot_selected_reads.R $sample
# # I have a problem with the dot_plot.pl on the server so locally:
# R version 4.1.2
# cd mnt/scratch/minion/
# gitHubDirectory="/home/ldelisle/mnt/home_scitas/softwares/scriptsForBoltEtAl2022/"
# path="$PWD/"
# samples=('542_nCATS' '320_nCATS')
# SLURM_ARRAY_TASK_ID=0
# sample=${samples[${SLURM_ARRAY_TASK_ID}]}
# pathResults=${path}/${sample}_results/
# cd $pathResults
# Rscript ${gitHubDirectory}/nCATS/scripts/dotplot_selected_reads.R $sample
