#!/bin/bash

#--- directories ---#
resultdir="/proj/yunligrp/users/chanhwa/pwas/whi_pwas_freeze9/result"
protein_list="/proj/yunligrp/users/chanhwa/pwas/whi_pwas_freeze9/data/proteindata/protein_list.txt"

#--- p-value thresholds used ---#
pval_thresh=(0.001 0.0005 0.0001 0.00005 0.00001 0.000005 0.000001)

#--- Generate header ---#
python -c 'import sys; header=["protein", "cisonly"] + ["proposed_" + pval for pval in sys.argv[1:]]; print("\t".join(header))' "${pval_thresh[@]}" > correlations.txt

#--- Read rsq values ---#
while read protein;do
  
  rsqs=(${protein})
  
  # Cis only method
  cisonly_rsq=$(cat   ${resultdir}/${protein}/cisonly/test.${protein}.log | grep 'Pearson R2 of predicted and observed protein level is' | grep -wo '0\.[0-9]*')

  [ -z $cisonly_rsq ] && cisonly_rsq="NA"
  
  rsqs+=(${cisonly_rsq})
  
  for pval in ${pval_thresh[@]};do  
  
    # Proposed method pval_thresh = 0.001
    proposed_rsq=$(cat   ${resultdir}/${protein}/proposed_pval_${pval}/test.${protein}.log | grep 'Pearson R2 of predicted and observed protein level is' | grep -wo '0\.[0-9]*')
  
    [ -z $proposed_rsq ] && proposed_rsq="NA"
  
  rsqs+=(${proposed_rsq})
  
  done
  
  python -c 'import sys; print("\t".join(sys.argv[1:]))' "${rsqs[@]}" >> correlations.txt

done < ${protein_list}
