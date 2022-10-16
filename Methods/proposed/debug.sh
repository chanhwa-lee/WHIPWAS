# Remove erroneous directories

resultdir="/proj/yunligrp/users/chanhwa/pwas/whi_pwas_freeze9/result"
protein_list="/proj/yunligrp/users/chanhwa/pwas/whi_pwas_freeze9/data/proteindata/protein_list.txt"
pval_thresh=0.01

while read protein;do

  outdir="${resultdir}/${protein}/proposed_pval_${pval_thresh}"

  rm -r ${outdir}
  
done < ${protein_list}