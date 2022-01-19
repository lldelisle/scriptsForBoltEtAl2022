# HiCExplorer version 3.6
pathForCool=/scratch/ldelisle/Bolt2022/cHi-C/
cd $pathForCool
for tis in PFL DFL; do
  #542-wt
  hicCompareMatrices \
    -m "542_E12.5_${tis}_CHIC_542mapping_5kb.cool" "wt_E12.5_${tis}_CHIC_542mapping_5kb.cool" \
    --operation diff \
    --outFileName "../542_minus_wt_E12.5_${tis}_CHIC_542mapping_5kb.cool"
done
