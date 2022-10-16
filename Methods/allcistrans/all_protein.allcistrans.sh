while read -r protein chr begin end; do

  echo Protein = $protein CHR = $chr Begin = $begin End = $end
  
  proteindir="/proj/yunligrp/users/chanhwa/pwas/whi_pwas_neat/result/${protein}"

  outdir="${proteindir}/allcistrans"
  
  rm -r ${outdir}

  mkdir -p ${outdir}
  
  sbatch -J ${protein}_allcistrans -t 300 --mem=16g --out ${outdir}/${protein}.allcistrans.log --wrap="bash allcistransmain.sh ${protein} $chr $begin $end"
  
done < "/proj/yunligrp/users/chanhwa/pwas/whi_pwas_neat/allcistrans/proteins_info.txt"