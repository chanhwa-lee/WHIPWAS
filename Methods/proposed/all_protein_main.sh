resultdir="/proj/yunligrp/users/chanhwa/pwas/whi_pwas_freeze9/result"
protein_list="/proj/yunligrp/users/chanhwa/pwas/whi_pwas_freeze9/data/proteindata/protein_list.txt"
pval_thresh=$1

while read protein;do

  outdir="${resultdir}/${protein}/proposed_pval_${pval_thresh}"

  mkdir -p ${outdir}
  
  sbatch -J ${protein}.pwas -o ${outdir}/${protein}.pwas.log -t 1- \
         --wrap="bash main.sh ${protein} ${pval_thresh}"
  
done < ${protein_list}