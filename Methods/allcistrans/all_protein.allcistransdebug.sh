while read -r protein chr begin end; do

  echo Protein = $protein CHR = $chr Begin = $begin End = $end
  
  proteindir="/proj/yunligrp/users/chanhwa/pwas/whi_pwas_neat/result/${protein}"

  outdir="${proteindir}/allcistrans"

  mkdir -p ${outdir}
  
  sbatch -J ${protein}_allcistrans -t 300 --mem=16g --out ${outdir}/${protein}.allcistransdebug.log --wrap="bash allcistransdebugmain.sh ${protein} $chr $begin $end"
  
done < "/proj/yunligrp/users/chanhwa/pwas/whi_pwas_neat/allcistrans/proteins_info.txt"