#!/bin/bash

#SBATCH -o slurm-%x-%A_%2a.out
#SBATCH -e slurm-%x-%A_%2a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lucille.delisle@epfl.ch
#SBATCH --nodes 1
#SBATCH --mem 150G
#SBATCH --cpus-per-task 36
#SBATCH --time 3-00:00:00
#SBATCH --array=0-1
#SBATCH --job-name basecall
#SBATCH --chdir /scratch/ldelisle/minion/

path="$PWD/"

samples=('542_nCATS' '320_nCATS')

sample=${samples[${SLURM_ARRAY_TASK_ID}]}

mkdir -p guppy_output_${sample}
~/softwares/ont-guppy-cpu/bin/guppy_basecaller --version
# : Guppy Basecalling Software, (C) Oxford Nanopore Technologies, Limited. Version 5.0.16+b9fcd7b
~/softwares/ont-guppy-cpu/bin/guppy_basecaller  -i  ${path}/${sample}/ \
  -s guppy_output_${sample}/ --flowcell FLO-MIN106 --kit SQK-LSK109 --num_callers 9 \
  --cpu_threads_per_caller 4 --fast5_out

mkdir -p basecalled_fast5_from_calling_${sample}
ln -s $PWD/guppy_output_${sample}*/workspace/*fast5 basecalled_fast5_from_calling_${sample}/

mkdir -p basecalled_fast5_single_${sample}
multi_to_single_fast5 -i basecalled_fast5_from_calling_${sample} -s basecalled_fast5_single_${sample}/ -t 36

mkdir -p basecalled_fast5_combined_${sample}
single_to_multi_fast5 -i basecalled_fast5_single_${sample}/ -s basecalled_fast5_combined_${sample}/ -t 36 -f ${sample} -n 1000000 --recursive
