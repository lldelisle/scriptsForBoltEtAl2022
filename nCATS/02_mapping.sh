#!/bin/bash

#SBATCH -o slurm-%x-%A_%2a.out
#SBATCH -e slurm-%x-%A_%2a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lucille.delisle@epfl.ch
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --mem 50G 
#SBATCH --cpus-per-task 10
#SBATCH --time 24:00:00
#SBATCH --array=0-1
#SBATCH --job-name nCATS_mapping
#SBATCH --chdir /scratch/ldelisle/minion/

gitHubDirectory=$1

samples=('542_nCATS' '320_nCATS')

path="$PWD/"
sample=${samples[${SLURM_ARRAY_TASK_ID}]}

nThreads=10

# First map on mm10:
genome="mm10"
pathForFasta="/home/ldelisle/genomes/fasta/${genome}.fa"
pathForMinimap2Index="/scratch/ldelisle/genomes/minimap2/${genome}.mmi"
pathForMinimap2IndexStorage="/work/updub/minimap2/${genome}.mmi"

mkdir -p $(dirname $pathForMinimap2Index)

module purge
module load gcc/8.4.0
module load samtools/1.10
module load bedtools2/2.27.1
module load python/3.7.7

pathResults=${path}/${sample}_results/

mkdir -p $pathResults
echo $sample
cd $pathResults

mkdir -p combinedFastq
fastq="${pathResults}/combinedFastq/${sample}.fastq.gz"
if [ -e $fastq ]; then
  echo "The fastq already exists"
else
  cat ${path}/guppy_output_${sample}*/pass/*.fastq | gzip > $fastq
fi

if [ ! -e $pathForMinimap2Index ]; then
  if [ -e $pathForMinimap2IndexStorage ]; then
    cp -r $pathForMinimap2IndexStorage $pathForMinimap2Index
  else
    ~/softwares/minimap2-2.15_x64-linux/minimap2 -t $nThreads -d $pathForMinimap2Index $pathForFasta
    cp -r $pathForMinimap2Index $pathForMinimap2IndexStorage
  fi
fi

if [ ! -e ${sample}_algn.sam ]; then
  ~/softwares/minimap2-2.15_x64-linux/minimap2 -t $nThreads -ax map-ont $pathForMinimap2Index $fastq > ${sample}_algn.sam
fi

if [ ! -e ${sample}_mappingStats_det.txt ]; then
  python ${gitHubDirectory}/nCATS/scripts/mappingStatsNanopore.py --output ${sample}_mappingStats.txt --outputDetail ${sample}_mappingStats_det.txt ${sample}_algn.sam
fi

# Due to a design flaw, BAM does not work with CIGAR strings with >65535 operations (SAM and CRAM work). 
# So this may not work
# Here it worked
if [ ! -e ${sample}_mapped_sorted.bam ]; then
  samtools sort --threads $nThreads -o ${sample}_mapped_sorted.bam ${sample}_algn.sam
fi

if [ ! -e ${pathForFasta}.fai ]; then
  samtools faidx $pathForFasta
fi

# We remove the non primary alignments
if [ ! -e  ${sample}_mapped_sorted_filtered2.bam.bai ]; then
  samtools view -b ${sample}_mapped_sorted.bam --threads $nThreads -F0x104 > ${sample}_mapped_sorted_filtered2.bam
  samtools index ${sample}_mapped_sorted_filtered2.bam
fi

if [ ! -e ${sample}_mapped_sorted_filtered2.bed.gz ]; then
  bedtools bamtobed -i ${sample}_mapped_sorted_filtered2.bam | gzip > ${sample}_mapped_sorted_filtered2.bed.gz
fi

# Region of the transgene
samtools view ${sample}_mapped_sorted_filtered2.bam chr2:74073413-74076528 | cut -f 1 > targetted_reads.txt

# Region of insertion
samtools view ${sample}_mapped_sorted_filtered2.bam chr2:75262998-75286118 | cut -f 1 >> targetted_reads.txt

cat targetted_reads.txt | sort | uniq > targetted_reads_uniq.txt

# seqtk version 1.3.0
/home/ldelisle/softwares/seqtk/seqtk subseq combinedFastq/${sample}.fastq.gz targetted_reads_uniq.txt > targetted_reads.fastq
/home/ldelisle/softwares/seqtk/seqtk seq -A targetted_reads.fastq > targetted_reads.fa

# Get the lengths:
cat targetted_reads.fa | awk '$0 ~ ">" {if(NR>1){print c}; c=0;printf substr($1,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > targetted_reads.sizes

targetted_reads=$(cat targetted_reads.txt | tr "\n" ",")
zcat ${sample}_mapped_sorted_filtered2.bed.gz | awk -v ir=$targetted_reads 'BEGIN{
  split(ir,myreads,",")
  for (i in myreads){
    reads[myreads[i]] = 1
  }}
  {if($4 in reads){print}}' | gzip > ${sample}_mapped_sorted_filtered2_targetted.bed.gz

mkdir -p ${gitHubDirectory}/nCATS/outputs/
awk '
BEGIN{
  print "read_name\tsize\tmapped_size\tstatus"
}
{
  if(NR == FNR){
    ## If this is the first file
    sizes[$1] = $2
  } else {
    # This is the second file
    mapped_size = $3 - $2
    if (mapped_size < 1.05 * sizes[$4] && mapped_size > 0.95 * sizes[$4]){
      print($4"\t"sizes[$4]"\t"mapped_size"\tfully_mapped")
    } else {
      print($4"\t"sizes[$4]"\t"mapped_size"\tTOCHECK")
    }
  }
}' targetted_reads.sizes <(zcat ${sample}_mapped_sorted_filtered2_targetted.bed.gz) > ${gitHubDirectory}/nCATS/outputs/${sample}_targetted_reads_status.txt

cp targetted_reads.fa ${gitHubDirectory}/nCATS/outputs/${sample}_targetted_reads.fa
