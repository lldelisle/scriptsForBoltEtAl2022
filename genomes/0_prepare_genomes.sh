#!/bin/bash

#SBATCH -o slurm-%x-%A.out
#SBATCH -e slurm-%x-%A.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lucille.delisle@epfl.ch
#SBATCH --nodes 1
#SBATCH --mem 80G
#SBATCH --cpus-per-task 1
#SBATCH --time 48:00:00
#SBATCH --job-name preparegenome
#SBATCH --chdir /home/ldelisle/genomes/fasta/

gitHubDirectory=$1

# Get mm10 from UCSC:
wget http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz
gunzip mm10.fa.gz

# To prepare the mutant genome creation
# Get the mutant chr2
wget "https://zenodo.org/record/5785457/files/chr2_II1TDOM-542.fa.gz?download=1" -O chr2_542.fa.gz
gunzip chr2_542.fa.gz
# Get all chrs from mm10 except chr2 using seqtk
cat mm10.fa | grep ">" | grep -v chr2 | sed 's/>//' > listOfChrs.txt
# I reformat them:
# This can be long:
seqtk seq -U mm10.fa | seqtk subseq -l 60 - listOfChrs.txt > allChrsIncludingContigsExceptchr2.fa
cat chr2_542.fa allChrsIncludingContigsExceptchr2.fa > mm10_II1TDOMlacZr4-CB542.fa

# We also extract the region around the transgene:
echo -e "chr2\t75260000\t75293000" > target.bed
seqtk subseq -l 60 mm10_II1TDOMlacZr4-CB542.fa target.bed > ${gitHubDirectory}/nCATS/fasta/542_around_Tg.fa
