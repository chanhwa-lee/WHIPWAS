ped_path="/proj/yunligrp/users/chanhwa/pwas/WHI_PWAS/Step_3_divide_samples"

assays=("cardiom" "cardio2" "cardio3" "inflam" "oncol" "neurol")

for assay in ${assays[@]};do
      
  protein_list=($(head ${ped_path}/WHI_olink_gwas.${assay}.ped -n 1 | cut -f28-))
  # NOTE: The number "28" of -f28- is from the fact that the number of non-protein fields(IDs, gender, age, array dummy, plate dummy) is 27.
  
  for protein in ${protein_list[@]};do

    reason=$(cat ${protein}/${protein}_cis_prediction.log | grep 'Less than 2 non-zero coefficients so no model reported' | grep -wo 'Less')
    echo ${#reason}
    
    if [ ${#reason} -eq 4 ]
    then
      echo For $protein, reason has length 4
      echo $protein >> no_corr_reason.txt
    fi
        
  done

done