# toolshed.g2.bx.psu.edu/repos/lparsons/cutadapt/cutadapt/1.16.8
# command_version:1.16

cutadapt  -j ${GALAXY_SLOTS:-1}     -a 'TruSeq R1'='AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'         -A 'TruSeq R2'='AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT'      --output='out1.fq.gz' --paired-output='out2.fq.gz'  --error-rate=0.1 --times=1 --overlap=3        --minimum-length=15 --pair-filter=any   --quality-cutoff=30   'CHIC_R1.fq.gz' 'CHIC_R2.fq.gz'  > report.txt

# toolshed.g2.bx.psu.edu/repos/bgruening/hicup_hicup/hicup_hicup/0.6.1.0
# command_version:HiCUP v0.6.1
BOWTIE_PATH_BASH="$(which bowtie2)" 
hicup_digester --re1 '^GATC' --genome 'mm10_II1TDOM' mm10_II1TDOMlacZr4-CB542.fasta.fasta 
mv *Digest_* digest_file.txt 
hicup --zip --threads ${GALAXY_SLOTS:-1} --digest digest_file.txt --index '/data/galaxy/galaxy/var/tool-data/mm10_II1TDOMlacZr4-CB542/bowtie2_index/mm10_II1TDOMlacZr4-CB542/mm10_II1TDOMlacZr4-CB542' --bowtie2 $BOWTIE_PATH_BASH --keep ou1.fq.gz out2.fq.gz  

# testtoolshed.g2.bx.psu.edu/repos/lldelisle/from_hicup_to_juicebox/fromHicupToJuicebox/0.0.2
# This tool can be downloaded from https://testtoolshed.g2.bx.psu.edu/repository/download?repository_id=be5040251cd4afb7&changeset_revision=44365a4feb3b&file_type=gz
python /data/galaxy/galaxy/var/shed_tools/testtoolshed.g2.bx.psu.edu/repos/lldelisle/from_hicup_to_juicebox/552eb7782435/from_hicup_to_juicebox/fromHicupToJuicebox.py --fragmentFile digest_file.txt --colForChr 1 --colForStart 2 --colForEnd 3 --colForID 4 --lineToSkipInFragmentFile 2 --methodForFrag hicup --useMid --output validPairs.txt "HiCUP data result.qname_sorted.bam" 
bash /data/galaxy/galaxy/var/shed_tools/testtoolshed.g2.bx.psu.edu/repos/lldelisle/from_hicup_to_juicebox/552eb7782435/from_hicup_to_juicebox/switchAndSort.sh validPairs.txt "validPairs file with midFrag positions.tabular"

# Filtering with c10>=30 and c11>=30 with filtering tool from galaxy on validPairs_sorted.txt

# Filtering with (c3=='chr2' and c4<77000000 and c4>72402000) and (c7=="chr2" and c8<77000000 and c8>72402000) with filtering tool from galaxy on the previous output.txt

# testtoolshed.g2.bx.psu.edu/repos/lldelisle/cooler/cooler_csort_pairix/0.0.1
# command_version:cooler, version 0.7.4

cooler csort -i tabix -c1 3 -c2 7 -p1 4 -p2 8 -o tabix_files.tabix "both pairs MAPQ30 and in captured region.tabular" mm10_II1TDOMlacZr4-CB542.sizes.gg

# testtoolshed.g2.bx.psu.edu/repos/lldelisle/cooler/cooler_makebins/0.0.1
# command_version:cooler, version 0.7.4

cooler makebins -o 5kb_bins.bed mm10_II1TDOMlacZr4-CB542.sizes.gg 5000

# testtoolshed.g2.bx.psu.edu/repos/lldelisle/cooler/cooler_cload_pairix/0.0.1
# command_version:cooler, version 0.7.4

cooler cload tabix --assembly mm10 -c2 7 -p2 8 5kb_bins.bed tabix_files.tabix raw_cool_file_5kb.cool

# testtoolshed.g2.bx.psu.edu/repos/lldelisle/cooler/cooler_balance/0.0.1
# command_version:cooler, version 0.7.4

cp raw_cool_file_5kb.cool balanced_cool_file_5kb.cool
cooler balance --mad-max 5 --min-nnz 10 --min-count 0 --ignore-diags 2 --tol 1e-05 --max-iters 200 --cis-only -f balanced_cool_file_5kb.cool
