genome is either mm10 or the mutant genome II1
# toolshed.g2.bx.psu.edu/repos/lparsons/cutadapt/cutadapt/1.16.8
# command_version:1.16
cutadapt  -j ${GALAXY_SLOTS:-1}     -a 'Truseq R1'='GATCGGAAGAGCACACGTCTGAACTCCAGTCAC'         -A 'TruSeq R2'='GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT'      --output='out1.fq.gz' --paired-output='out2.fq.gz'  --error-rate=0.1 --times=1 --overlap=3        --minimum-length=15 --pair-filter=any   --quality-cutoff=30   'CUTandRUN_R1.fq.gz' 'CUTandRUN_R2.fq.gz'  > report.txt

# toolshed.g2.bx.psu.edu/repos/devteam/bowtie2/bowtie2/2.4.2+galaxy0
# command_version:2.4.1
bowtie2  -p ${GALAXY_SLOTS:-4}  -x 'genome'   -1 'input_f.fastq.gz' -2 'input_r.fastq.gz' -I 0 -X 1000 --fr --no-mixed --no-discordant --dovetail                --very-sensitive   2> 'mapping stats.txt'  | samtools sort --no-PG -@${GALAXY_SLOTS:-2} -T "${TMPDIR:-.}" -O bam -o 'bowtie2 output (BAM).bam'

# toolshed.g2.bx.psu.edu/repos/devteam/samtool_filter2/samtool_filter2/1.8+galaxy1
# command_version:
samtools view -o 'filtered BAM.bam' -h   -b  -q 30 -f 0x2 input.bam

# toolshed.g2.bx.psu.edu/repos/devteam/picard/picard_MarkDuplicates/2.18.2.2
# command_version:
_JAVA_OPTIONS=${_JAVA_OPTIONS:-'-Xmx2048m -Xms256m'} 
export _JAVA_OPTIONS 
picard MarkDuplicates  INPUT='filtered BAM.bam' OUTPUT='BAM filtered rmDup.bam'  METRICS_FILE='MarkDuplicates metrics.txt'  REMOVE_DUPLICATES='false' ASSUME_SORTED='true'  DUPLICATE_SCORING_STRATEGY='SUM_OF_BASE_QUALITIES'  OPTICAL_DUPLICATE_PIXEL_DISTANCE='100'   VALIDATION_STRINGENCY='LENIENT' QUIET=true VERBOSITY=ERROR  `if [ -n "$TMPDIR" ] ; then echo 'TMP_DIR=$TMPDIR' ; else if [ -n "$TEMP" ] ; then echo 'TMP_DIR=$TEMP' ; fi ; fi`

# toolshed.g2.bx.psu.edu/repos/iuc/bedtools/bedtools_bamtobed/2.29.2
# command_version:bedtools v2.29.2
bedtools bamtobed    -i 'BAM filtered rmDup.bam' > 'BED filtered rmDup.bed'

# toolshed.g2.bx.psu.edu/repos/iuc/macs2/macs2_callpeak/2.1.1.20160309.6
# command_version:macs2 2.1.1.20160309
macs2 callpeak   -t 'BED filtered rmDup.bed'  --name CUTandRUN    --format BED   --gsize '1870000000'           --call-summits  --keep-dup 'all'   --bdg  --qvalue '0.05'  --nomodel --extsize '200' --shift '-100'  2>&1 > macs2_stderr

# wig_to_bigWig
# command_version:
grep -v "^track" 'MACS2 treatment coverage.bedgraph' | wigToBigWig stdin 'genome.len' 'Coverage from MACS2 (bigwig).bigwig' -clip 2>&1 || echo "Error running wigToBigWig." >&2
