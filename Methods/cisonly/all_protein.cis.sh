while read -r protein chr begin end; do

  echo Protein = $protein CHR = $chr Begin = $begin End = $end

  outdir="/proj/yunligrp/users/chanhwa/pwas/whi_pwas_freeze9/result/${protein}/cisonly"

  mkdir -p ${outdir}
  
  sbatch -J ${protein}.cis -t 45 --mem=50g --out ${outdir}/${protein}.cisonly.log --wrap="bash main.sh ${protein} $chr $begin $end"
  
done < proteins_info.txt