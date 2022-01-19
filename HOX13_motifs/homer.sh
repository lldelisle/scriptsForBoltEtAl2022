# Homer version 4.11
# Galaxy Tool ID: 	toolshed.g2.bx.psu.edu/repos/iuc/homer_findmotifsgenome/homer_findMotifsGenome/4.11+galaxy2
findMotifsGenome.pl 'ChIPorCUTandRUN.narrowPeak' 'mm10_UCSC.fa' 'ChIPorCUTandRUN_motif' \
    -preparsedDir '/data/galaxy/galaxy/var/tool-data/mm10_UCSC/homer_preparse/mm10_UCSC_50' \
    -size '50' -len '8,10,12' -S 25 -mis 2    -mset auto \
    -gc  -local 0 -redundant 2.0 -maxN 0.7   -nlen '3' -nmax '160'  -e '0.0'  -minlp '-10.0'

# CB_II1 contains chr2:74075312-74075843
# Galaxy Tool ID: 	toolshed.g2.bx.psu.edu/repos/iuc/homer_scanmotifgenomewide/homer_scanMotifGenomeWide/4.11+galaxy2
scanMotifGenomeWide.pl 'ChIPorCUTandRUN_motif1.motif' 'CB_II1.fa' -bed > 'ChIPorCUTandRUN_motif1_II1.bed'
# To go back to mm10 and 0-based:
awk '{$1 = "chr2"; $2 += 74075311 - 1; $3 += 74075311; print}' 'ChIPorCUTandRUN_motif1_II1.bed' > 'ChIPorCUTandRUN_motif1.bed'

# To get the best motif:
cat HOX13_motifs/outputs/*motif1.bed | awk -v OFS="\t" '
{
  cur_chr = $1
  cur_start = $2
  cur_end = $3
  cur_mid = sprintf("%d", (cur_end + cur_start) / 2)
  cur_motif = $4
  cur_score = $5
  cur_strand = $6
  if (cur_mid in scores_per_region){
    if (scores_per_region[cur_mid] < cur_score){
      scores_per_region[cur_mid] = cur_score
      infos_per_region[cur_mid]["chr"] = cur_chr
      infos_per_region[cur_mid]["start"] = cur_start
      infos_per_region[cur_mid]["end"] = cur_end
      infos_per_region[cur_mid]["motif"] = cur_motif
      infos_per_region[cur_mid]["score"] = cur_score
      infos_per_region[cur_mid]["strand"] = cur_strand
    }
  } else {
    scores_per_region[cur_mid] = cur_score
    infos_per_region[cur_mid]["chr"] = cur_chr
    infos_per_region[cur_mid]["start"] = cur_start
    infos_per_region[cur_mid]["end"] = cur_end
    infos_per_region[cur_mid]["motif"] = cur_motif
    infos_per_region[cur_mid]["score"] = cur_score
    infos_per_region[cur_mid]["strand"] = cur_strand
  }
}
END{
  for (cur_mid in infos_per_region){
    printf "%s\t%d\t%d\t%s\t%f\t%s\n",
    infos_per_region[cur_mid]["chr"], \
    infos_per_region[cur_mid]["start"], \
    infos_per_region[cur_mid]["end"], \
    infos_per_region[cur_mid]["motif"], \
    infos_per_region[cur_mid]["score"], \
    infos_per_region[cur_mid]["strand"]
  }
}' > HOX13_motifs/outputs/best.bed

# The three motifs within the highlighted region in Figure 2A:
# First is TTTATDGG (SRR3498935_wt_E11_FL_HOXD13)
# Second is SCTRTAAA (wt_E12.5_DFL_HOXA13)
# Last is TTTATRGS (wt_E12.5_DFL_HOXD13_rep3)
