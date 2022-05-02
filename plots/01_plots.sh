# Using pygenometracks version Bolt2022

gitHubDirectory="/home/ldelisle/softwares/scriptsForBoltEtAl2022/"
pathWithCool="/scratch/ldelisle/Bolt2022/cHi-C/"
pathWithCR="/scratch/ldelisle/Bolt2022/CUTandRUN/"
pathWithATAC="/scratch/ldelisle/Bolt2022/ATAC/"
pathWithChIP="/scratch/ldelisle/Bolt2022/ChIP/"
plottingDirectory="/scratch/ldelisle/Bolt2022/plots/"
mkdir -p $plottingDirectory

# Figure 1A:
ini_file=${plottingDirectory}/fig1a.ini
echo "[scalebar]
file_type = scalebar
line_width = 1.5
size = 250000
where = top
scalebar_start_position = 73950000

[ATAC wt DFL Bolt2021 mm10Genome]
file = ${pathWithATAC}/wt_E12.5_DFL_ATAC_Bolt2021_mm10mapping.bw
title = ATAC wt DFL
height = 3
color = #999999
min_value=0
max_value=100

[spacer]

[wt_E12.5_DFL_HOXA13_rep2]
file = ${pathWithCR}/wt_E12.5_DFL_HOXA13_rep2_mm10mapping.bw
title = HOXA13
height = 3
color = #56b4e9ff
min_value=0
max_value=100

[spacer]

[wt_E12.5_DFL_HOXD13_rep2]
file = ${pathWithCR}/wt_E12.5_DFL_HOXD13_rep2_mm10mapping.bw
title = HOXD13
height = 3
color = #0072b2ff
min_value=0
max_value=100

[spacer]

[features]
file = ${gitHubDirectory}/annotations/HoxD_Elements_mm10_colored.bed
title = features
display = collapsed
color = bed_rgb
line_width = 2
border_color = bed_rgb
labels = false

[spacer]

[features labels]
file = ${gitHubDirectory}/annotations/HoxD_Elements_mm10_start.bed
display = collapsed
color = none
line_width = 0
labels = true
" > ${ini_file}

pgt --tracks ${ini_file} \
    -o ${ini_file/.ini/.pdf} --region chr2:73950000-75655000 \
    --dpi 300 --fontSize 12 --trackLabelFraction 0.1


# Figure 1B:
# Extract summits:
zcat ${pathWithCR}/wt_E12.5_DFL_HOXA13_rep2_mm10mapping.narrowPeak.gz | awk -v OFS="\t" '{print $1, $2 + $10, $2 + $10 + 1}' > ${plottingDirectory}/wt_E12.5_DFL_HOXA13_rep2_mm10mapping_summits.bed
zcat ${pathWithCR}/wt_E12.5_DFL_HOXD13_rep2_mm10mapping.narrowPeak.gz | awk -v OFS="\t" '{print $1, $2 + $10, $2 + $10 + 1}' > ${plottingDirectory}/wt_E12.5_DFL_HOXD13_rep2_mm10mapping_summits.bed
ini_file=${plottingDirectory}/fig1b.ini
echo "[scalebar]
file_type = scalebar
height = 0.5
line_width = 2
where = bottom
size = 500
scalebar_start_position = 74072912

[spacer]

[H3K27ac wt DFL ERC2017]
file = ${pathWithChIP}/SRR5855214_wt_E12_DFL_H3K27ac.bw
title = H3K27ac wt DFL
height = 3
color = #999999
min_value = 0
max_value = 200

[spacer]

[ATAC wt DFL Bolt2021 mm10Genome]
file = ${pathWithATAC}/wt_E12.5_DFL_ATAC_Bolt2021_mm10mapping.bw
title = ATAC wt DFL
height = 3
color = #999999
min_value = 0
max_value = 300

[spacer]

[wt_E12.5_DFL_HOXA13_rep2_mm10mapping]
file = ${pathWithCR}/wt_E12.5_DFL_HOXA13_rep2_mm10mapping.bw
title = HOXA13
height = 3
color = #56b4e9ff
min_value = 0
max_value = 300

[HOXA13 Sheth]
file = ${pathWithChIP}/SRR3498934_wt_E11_FL_HOXA13.bw
title = HOXA13 Sheth
color = #8b8b8bff
type = line
line_width = 1
overlay_previous = share-y

[summit]
file = ${plottingDirectory}/wt_E12.5_DFL_HOXA13_rep2_mm10mapping_summits.bed
title = summit
display = stacked
color = #56b4e9ff
line_width = 4
border_color = #56b4e9ff
labels = false

[spacer]

[wt_E12.5_DFL_HOXD13_rep2]
file = ${pathWithCR}/wt_E12.5_DFL_HOXD13_rep2_mm10mapping.bw
title = HOXD13
height = 3
color = #0072b2ff
min_value = 0
max_value = 300

[HOXD13 Sheth]
file = ${pathWithChIP}/SRR3498935_wt_E11_FL_HOXD13.bw
title = HOXD13 Sheth
color = #8b8b8bff
type = line
line_width = 1
overlay_previous = share-y

[summit]
file = ${plottingDirectory}/wt_E12.5_DFL_HOXD13_rep2_mm10mapping_summits.bed
title = summit
display = stacked
color = #0072b2ff
line_width = 4
border_color = #0072b2ff
labels = false

[spacer]

[features]
file = ${gitHubDirectory}/annotations/HoxD_Elements_mm10_colored.bed
labels_in_margin = true
display = stacked
height = 3
line_width = 0
color = bed_rgb

[highlight]
file = ${gitHubDirectory}/annotations/aroundHox13Motifs.bed
type = vhighlight
color = #dddadaff
" > ${ini_file}

pgt --tracks ${ini_file} \
  -o ${ini_file/.ini/.pdf} --region chr2:74072912-74077028 \
  --dpi 300 --fontSize 12 --trackLabelFraction 0.1 --plotWidth 15

# Figure 1SA:
ini_file=${plottingDirectory}/fig1sa.ini
echo "[scalebar]
file_type = scalebar
height = 0.5
line_width = 2
size = 1000

[wt_E12.5_DFL_HOXA13_rep2_mm10mapping]
file = ${pathWithCR}/wt_E12.5_DFL_HOXA13_rep2_mm10mapping.bw
title = HOXA13
height = 3
color = #56b4e9ff
min_value = 0

[HOXA13 Sheth]
file = ${pathWithChIP}/SRR3498934_wt_E11_FL_HOXA13.bw
title = HOXA13 Sheth
color = #8b8b8bff
type = line
line_width = 1
min_value = 0
overlay_previous = yes

[spacer]

[wt_E12.5_DFL_HOXD13_rep2]
file = ${pathWithCR}/wt_E12.5_DFL_HOXD13_rep2_mm10mapping.bw
title = HOXD13
height = 3
color = #0072b2ff
min_value = 0

[HOXD13 Sheth]
file = ${pathWithChIP}/SRR3498935_wt_E11_FL_HOXD13.bw
title = HOXD13 Sheth
color = #8b8b8bff
type = line
line_width = 1
min_value = 0
overlay_previous = yes

[spacer]

[features]
file = ${gitHubDirectory}/annotations/HoxD_Elements_mm10_colored.bed
labels_in_margin = true
display = stacked
height = 3
line_width = 0
color = bed_rgb
" > ${ini_file}

pgt --tracks ${ini_file} \
  -o ${ini_file/.ini/.pdf} --region chr2:74072912-74077028 \
  --dpi 300 --fontSize 12 --trackLabelFraction 0.1 --plotWidth 15


# Figure 3Aleft:
ini_file=${plottingDirectory}/fig3aleft.ini
echo "[scalebar]
file_type = scalebar
height = 0.5
line_width = 2
size = 500
where = bottom
scalebar_start_position = 74074674
" > ${ini_file}
for tis in FB PFL; do
    echo "[ATAC wt ${tis} rep1]
file = ${pathWithATAC}/wt_E12.5_${tis}_ATAC_rep1_mm10mapping.bw
title = ATAC ${tis} rep1
height = 3
color = #999999
min_value = 0
max_value = 100

[ATAC wt ${tis} rep2]
file = ${pathWithATAC}/wt_E12.5_${tis}_ATAC_rep2_mm10mapping.bw
title = ATAC ${tis} rep2
color = black
type = line
line_width = 1
overlay_previous = share-y

[spacer]
height = 1
" >> ${ini_file}
done
tis="DFL"
echo "[ATAC wt ${tis} rep1]
file = ${pathWithATAC}/wt_E12.5_${tis}_ATAC_rep1_mm10mapping.bw
title = ATAC ${tis} rep1
height = 3
color = grey
min_value = 0
max_value = 100

[ATAC wt ${tis} rep2]
file = ${pathWithATAC}/wt_E12.5_${tis}_ATAC_rep2_mm10mapping.bw
title = ATAC ${tis} rep2
color = black
type = line
line_width = 1
min_value = 0
overlay_previous = yes

[spacer]
height = 1
" >> ${ini_file}

echo "[wt_E12.5_DFL_HOXA13_rep2_mm10mapping]
file = ${pathWithCR}/wt_E12.5_DFL_HOXA13_rep2_mm10mapping.bw
title = HOXA13 wt rep2
height = 3
color = #56b4e9ff
min_value = 0

[wt_E12.5_DFL_HOXA13_rep1_mm10mapping]
file = ${pathWithCR}/wt_E12.5_DFL_HOXA13_rep1_mm10mapping.bw
title = HOXA13 wt rep1
color = black
min_value = 0
type = line
line_width = 1
overlay_previous = yes

[spacer]
height = 1

[wt_E12.5_DFL_HOXD13_rep2_mm10mapping]
file = ${pathWithCR}/wt_E12.5_DFL_HOXD13_rep2_mm10mapping.bw
title = HOXD13 wt rep2
height = 3
color = #0072b2ff
min_value = 0
max_value = 300

[wt_E12.5_DFL_HOXD13_rep1_mm10mapping]
file = ${pathWithCR}/wt_E12.5_DFL_HOXD13_rep1_mm10mapping.bw
title = HOXD13 wt rep1
color = black
min_value = 0
type = line
line_width = 1
overlay_previous = yes

[spacer]

[features]
file = ${gitHubDirectory}/annotations/HoxD_Elements_mm10_colored.bed
labels_in_margin = true
display = stacked
height = 3
line_width = 0
color = bed_rgb

[HOX]
file = ${gitHubDirectory}/annotations/Hox13_motifs_mm10_colored_homer.bed
title = motifs
display = collapsed
line_width = 0.5
border_color = bed_rgb
color = bed_rgb
labels = false

[highlight]
file = ${gitHubDirectory}/annotations/aroundHox13Motifs.bed
type = vhighlight
color = #dddadaff
" >> ${ini_file}

pgt --tracks ${ini_file} \
  -o ${ini_file/.ini/.pdf} --region chr2:74074674-74076672 \
  --dpi 300 --fontSize 12 --trackLabelFraction 0.2 --plotWidth 7


# Figure 3Aright:
ini_file=${plottingDirectory}/fig3aright.ini
echo "[scalebar]
file_type = scalebar
height = 0.5
line_width = 2
size = 500
where = bottom
scalebar_start_position = 75268925
" > ${ini_file}
for tis in FB PFL DFL; do
    echo "[ATAC 542 ${tis} rep1]
file = ${pathWithATAC}/542_E12.5_${tis}_ATAC_rep1_542mapping.bw
title = ATAC ${tis} rep1
height = 3
color = #999999
min_value = 0
max_value = 100

[ATAC 542 ${tis} rep2]
file = ${pathWithATAC}/542_E12.5_${tis}_ATAC_rep2_542mapping.bw
title = ATAC ${tis} rep2
color = black
type = line
line_width = 1
overlay_previous = share-y

[spacer]
height = 1
" >> ${ini_file}
done

echo "[542_E12.5_DFL_HOXA13_rep2_542mapping]
file = ${pathWithCR}/542_E12.5_DFL_HOXA13_rep2_542mapping.bw
title = HOXA13 542 rep2
height = 3
color = #56b4e9ff
min_value = 0

[542_E12.5_DFL_HOXA13_rep1_542mapping]
file = ${pathWithCR}/542_E12.5_DFL_HOXA13_rep1_542mapping.bw
title = HOXA13 542 rep1
color = black
min_value = 0
type = line
line_width = 1
overlay_previous = yes

[spacer]
height = 1

[542_E11.5_DFL_HOXD13_542mapping]
file = ${pathWithCR}/542_E11.5_DFL_HOXD13_542mapping.bw
title = HOXD13 CB542 E11
height = 3
color = #0072b2ff
min_value = 0

[spacer]

[features]
file = ${gitHubDirectory}/annotations/CB542_annotations_colored.bed
title = features
display = collapsed
line_width = 0
color = bed_rgb
labels = false

[spacer]

[features]
file = ${gitHubDirectory}/annotations/CB542_annotations_colored.bed
labels_in_margin = true
display = stacked
line_width = 0
color = none
labels = true
height = 1.5

[spacer]

[HOX]
file = ${gitHubDirectory}/annotations/Hox13_motifs_542_colored_homer.bed
title = motifs
display = collapsed
line_width = 0.5
border_color = bed_rgb
color = bed_rgb
labels = false

[highlight]
file = ${gitHubDirectory}/annotations/aroundHox13Motifs_542.bed
type = vhighlight
color = #dddadaff
" >> ${ini_file}

pgt --tracks ${ini_file} \
  -o ${ini_file/.ini/.pdf} --region chr2:75268925-75270923 \
  --dpi 300 --fontSize 12 --trackLabelFraction 0.2 --plotWidth 7

# Figure 3B
ini_file=${plottingDirectory}/fig3b.ini
echo "
[features]
file = ${gitHubDirectory}/annotations/II1_mm10_colored.bed
labels_in_margin = true
display = stacked
line_width = 0
color = bed_rgb

[HOX]
file = ${gitHubDirectory}/annotations/Hox13_motifs_mm10_colored_homer.bed
display = collapsed
line_width = 0
color = bed_rgb
labels = false
overlay_previous = yes
" > ${ini_file}

pgt --tracks ${ini_file} \
  -o ${ini_file/.ini/.pdf} --region chr2:74075300-74075900 \
  --dpi 300 --fontSize 12 --trackLabelFraction 0.1 --plotWidth 15

# Figure 4A
wget http://ftp.ensembl.org/pub/release-102/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.chromosome.2.fa.gz
gunzip Mus_musculus.GRCm38.dna.chromosome.2.fa.gz
ini_file=${plottingDirectory}/fig4a.ini
echo "[homer motifs]
file = ${gitHubDirectory}/annotations/Hox13_motifs_mm10_colored_homer.bed
color = bed_rgb
labels = false
display = collapsed

[fasta]
file = $PWD/Mus_musculus.GRCm38.dna.chromosome.2.fa

[crispr guides]
file = ${gitHubDirectory}/annotations/CRISPR_guides_for_II1_TFBS.bed
display = collapsed
labels = false

[crispr guides labels]
file = ${gitHubDirectory}/annotations/CRISPR_guides_for_II1_TFBS_labels.bed
display = collapsed
overlay_previous = yes
line_width = 0
color = none

[vlines]
file = ${gitHubDirectory}/annotations/CRISPR_guides_for_II1_TFBS_labels.bed
type = vlines
" > ${ini_file}

pgt --tracks ${ini_file} \
  -o ${ini_file/.ini/.pdf} --region chr2:74,075,687-74,075,808 \
  --dpi 300 --fontSize 12 --trackLabelFraction 0.1 --plotWidth 15


# Figure 6A
ini_file=${plottingDirectory}/fig6a.ini
echo "" > ${ini_file}
for tis in PFL DFL; do
    echo "[542-wt ${tis}]
file = ${pathWithCool}/../542_minus_wt_E12.5_${tis}_CHIC_542mapping_5kb.cool
title = 542-wt ${tis} 5kb
depth = 1300000
show_masked_bins = false
colormap = bwr
transform = no
min_value = -4e-06
max_value = 4e-06

[spacer]
" >> ${ini_file}
done
echo "[annot]
file = ${gitHubDirectory}/annotations/HoxD_Elements_542_colored.bed
title = 542 HoxD elements
display = collapsed
color = bed_rgb
line_width = 2
border_color = bed_rgb
labels = false

[spacer]

[features labels]
file = ${gitHubDirectory}/annotations/HoxD_Elements_542_start.bed
display = collapsed
color = none
line_width = 0
labels = true
" >> ${ini_file}

pgt --tracks ${ini_file} \
  -o ${ini_file/.ini/.pdf} --region chr2:73950000-75655000 \
  --dpi 300 --fontSize 12 --trackLabelFraction 0.1

# Figure 6B
ini_file=${plottingDirectory}/fig6b.ini
tis="PFL"
echo "
[542 ${tis}]
file = ${pathWithCool}/542_E12.5_${tis}_CHIC_542mapping_5kb.cool
region2 = chr2:75255389-75285908
title = 542 ${tis} 5kb
show_masked_bins = false
colormap = YlGnBu
transform = no
min_value = 0.00
max_value = 0.01
file_type = hic_matrix_square

[Hoxd8_Tg]
file = ${gitHubDirectory}/annotations/Hoxd8_Tg.bedpe
links_type = squares
line_width = 0.5
color = red
overlay_previous = share-y
region2 = chr2:75255389-75285908

[spacer]

[scalebar]
file_type = scalebar
height = 0.5
line_width = 2
size = 10000
where = bottom
scalebar_start_position = 74638000

[ERC wt ${tis} H3K27ac]
file = ${pathWithChIP}/SRR5855215_wt_E12_${tis}_H3K27ac.bw
title = ERC wt ${tis} Ac
height = 3
color = grey
min_value = 0
max_value = 400

[spacer]

[Beccari wt ${tis} H3K27me3]
file = ${pathWithChIP}/SRR3168462_wt_E12_${tis}_H3K27me3.bw
title = Beccari wt ${tis} me3
height = 3
color = grey
min_value=0
max_value = 200

[spacer]

[ERC wt ${tis} CTCF]
file = ${pathWithChIP}/SRR5855221_wt_E12_${tis}_CTCF.bw
title = ERC wt ${tis} CTCF
height = 3
color = black
min_value = 0
max_value = 200

[spacer]

[annot]
file = ${gitHubDirectory}/annotations/HoxD_Elements_mm10_colored.bed
title = wt HoxD elements
display = collapsed
color = bed_rgb
line_width = 0
labels = false

[spacer]

[features labels]
file = ${gitHubDirectory}/annotations/HoxD_genes_mm10_start.bed
display = collapsed
fontstyle = italic
color = none
line_width = 0
labels = true
" > ${ini_file}

pgt --tracks ${ini_file} \
  -o ${ini_file/.ini/.pdf} --region chr2:74638000-74783000 \
  --dpi 300 --fontSize 12 --trackLabelFraction 0.2

# Figure 6C
ini_file=${plottingDirectory}/fig6c.ini
tis="DFL"
echo "[542 ${tis}]
file = ${pathWithCool}/542_E12.5_${tis}_CHIC_542mapping_5kb.cool
region2 = chr2:75255389-75285908
title = 542 ${tis} 5kb
show_masked_bins = false
colormap = YlGnBu
transform = no
min_value = 0.00
max_value = 0.01
file_type = hic_matrix_square

[Hoxd8_Tg]
file = ${gitHubDirectory}/annotations/Hoxd8_Tg.bedpe
links_type = squares
line_width = 0.5
color = red
overlay_previous = share-y
region2 = chr2:75255389-75285908

[spacer]

[scalebar]
file_type = scalebar
height = 0.5
line_width = 2
size = 10000
where = bottom
scalebar_start_position = 74638000

[spacer]

[ERC wt ${tis} H3K27ac]
file = ${pathWithChIP}/SRR5855214_wt_E12_${tis}_H3K27ac.bw
title = ERC wt ${tis} Ac
height = 3
color = grey
min_value = 0
max_value = 400

[spacer]

[Beccari wt ${tis} H3K27me3]
file = ${pathWithChIP}/SRR3168464_wt_E12_${tis}_H3K27me3.bw
title = Beccari wt ${tis} me3
height = 3
color = grey
min_value=0
max_value = 200

[spacer]

[ERC wt ${tis} CTCF]
file = ${pathWithChIP}/SRR5855220_wt_E12_${tis}_CTCF.bw
title = ERC wt ${tis} CTCF
height = 3
color = black
min_value = 0
max_value = 200

[spacer]

[annot]
file = ${gitHubDirectory}/annotations/HoxD_Elements_mm10_colored.bed
title = wt HoxD elements
display = collapsed
color = bed_rgb
line_width = 0
labels = false

[spacer]

[features labels]
file = ${gitHubDirectory}/annotations/HoxD_genes_mm10_start.bed
display = collapsed
fontstyle = italic
color = none
line_width = 0
labels = true
" > ${ini_file}

pgt --tracks ${ini_file} \
  -o ${ini_file/.ini/.pdf} --region chr2:74638000-74783000 \
  --dpi 300 --fontSize 12 --trackLabelFraction 0.2

# Figure 6S
ini_file=${plottingDirectory}/fig6s.ini
echo "" > ${ini_file}
for tis in PFL DFL; do 
    for geno in wt 542; do
        echo "[hic matrix]
file = ${pathWithCool}/${geno}_E12.5_${tis}_CHIC_542mapping_5kb.cool 
title = ${geno} ${tis} 5kb 542_mapping
depth = 1300000
show_masked_bins = false
colormap = YlGnBu
transform = no
min_value = 0.00
max_value = 0.01

[spacer]
" >> ${ini_file}
    done
done
echo "[annot]
file = ${gitHubDirectory}/annotations/HoxD_Elements_542_colored.bed
title = 542 HoxD elements
display = collapsed
color = bed_rgb
line_width = 2
border_color = bed_rgb
labels = false

[spacer]

[features labels]
file = ${gitHubDirectory}/annotations/HoxD_Elements_542_start.bed
display = collapsed
color = none
line_width = 0
labels = true
" >> ${ini_file}


pgt --tracks ${ini_file} \
  -o ${ini_file/.ini/.pdf} --region chr2:73950000-75655000 \
  --dpi 300 --fontSize 12 --trackLabelFraction 0.1

