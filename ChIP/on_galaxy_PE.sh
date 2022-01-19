# toolshed.g2.bx.psu.edu/repos/lparsons/cutadapt/cutadapt/1.16.8
# command_version:1.16
cutadapt  -j ${GALAXY_SLOTS:-1}     -a 'TrueSeq'='GATCGGAAGAGCACACGTCTGAACTCCAGTCAC'         -A 'TruSeq'='GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT'      --output='out1.fq.gz' --paired-output='out2.fq.gz'  --error-rate=0.1 --times=1 --overlap=3        --minimum-length=15 --pair-filter=any   --quality-cutoff=30   'ChIP_R1.fq.gz' 'ChIP_R2.fq.gz'  > report.txt

# toolshed.g2.bx.psu.edu/repos/devteam/bowtie2/bowtie2/2.4.2+galaxy0
# command_version:2.4.1
bowtie2  -p ${GALAXY_SLOTS:-4}  -x '/data/galaxy/galaxy/var/tool-data/mm10_UCSC/bowtie2_index/mm10_UCSC/mm10_UCSC'   -1 'out1.fq.gz' -2 'out2.fq.gz'                2> 'mapping stats.txt'  | samtools sort --no-PG -@${GALAXY_SLOTS:-2} -T "${TMPDIR:-.}" -O bam -o 'bowtie2 output (BAM).bam'

# toolshed.g2.bx.psu.edu/repos/devteam/samtool_filter2/samtool_filter2/1.8+galaxy1
# command_version:
samtools view -o 'filtered BAM.bam' -h   -b  -q 30 -f 0x2 'bowtie2 output (BAM).bam'

# toolshed.g2.bx.psu.edu/repos/iuc/macs2/macs2_callpeak/2.1.1.20160309.6
# command_version:macs2 2.1.1.20160309
macs2 callpeak   -t 'filtered BAM.bam'  --name SRR3498934    --format BAMPE   --gsize '1870000000'           --call-summits  --keep-dup '1'   --bdg  --qvalue '0.05'  --mfold '5' '50'  --bw '300'  2>&1 > macs2_stderr

# wig_to_bigWig
# command_version:
grep -v "^track" 'MACS2 treatment coverage.bedgraph' | wigToBigWig stdin '/data/galaxy/galaxy/var/tool-data/mm10_UCSC/len/mm10_UCSC.len' 'coverage from MACS2 (bigwig).bigwig' -clip 2>&1 || echo "Error running wigToBigWig." >&2
