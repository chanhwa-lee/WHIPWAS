#!/bin/bash

#----- Variables for cisonly method -----#

starttime=`date +%s`

protein=$1
chr=$2
begin=$3
end=$4

genodir="/proj/yunligrp/users/chanhwa/pwas/whi_pwas_freeze9/data/genotypedata"

peddir="/proj/yunligrp/users/chanhwa/pwas/whi_pwas_freeze9/data/proteindata"

tooldir="/proj/yunligrp/users/chanhwa/pwas/whi_pwas_freeze9/method/cisonly"

outdir="/proj/yunligrp/users/chanhwa/pwas/whi_pwas_freeze9/result/${protein}/cisonly"

#----- Makde output directory -----#

mkdir -p ${outdir}

cd ${outdir}


echo -e "----------------------------------------------------------------------------------------------"
echo -e "----------------------------------------------------------------------------------------------\n"
echo -e "        Protein ${protein} Cis variants only prediction EN model training and testing         \n" 
echo -e "                              `date`                                                          \n"
echo -e "----------------------------------------------------------------------------------------------"
echo -e "----------------------------------------------------------------------------------------------\n"

#------ Load Modules ------#

module add python
module add samtools
module add r/4.0.3


echo -e "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"
echo -e "                                     Training Step                                          \n" 
echo -e "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"

# Check whether Training data is already generated. If not, generate it.

if [ -f "${outdir}/train.${protein}.dose.txt" ]
then
  echo -e "--- Training data is already generated. Proceed to next step. ---\n" 

else


  #------ Training dose file ------#
    
  start=$(( $begin - 1000000 ))
  if [[ $start -lt 1 ]]; then  
    start=1
  fi
  finish=$(( $end + 1000000 ))  
  
  if(( ${chr}==23 )); then chr=X; fi
      
  # Subset only cis variants from genotype data file
  bcftools view \
    -t ${chr}:$start-$finish \
    ${genodir}/train/WHI.olink.chr${chr}.dose.nochr.train.vcf.gz \
    -Oz -o ${outdir}/train.filtered.vcf.gz 

  # Generate header
  echo "CHROM POS REF ALT ID" $(bcftools query -l "${outdir}/train.filtered.vcf.gz") | sed 's/ /\t/g' > ${outdir}/header.txt

  # Extract Hard Genotype Calling (GT field) of cis variants from genotype data file
  bcftools query \
    -f '%CHROM\t%POS\t%REF\t%ALT{0}\t%ID[\t%GT]\n' \
    "${outdir}/train.filtered.vcf.gz" > ${outdir}/chr${chr}.gt.txt
  
  # Transform GT field (0/0, 0/1, 1/1) into dose (0, 1, 2)
  sed 's/0\/0/0/g' ${outdir}/chr${chr}.gt.txt | sed 's/0\/1/1/g' | sed 's/1\/1/2/g' > ${outdir}/chr${chr}.dose.txt
  
  echo -e "--- Training dose of chromosome ${chr} file generated ---\n"

  # Paste header and training dose file over chromosomes
  cat ${outdir}/header.txt ${outdir}/chr*.dose.txt > ${outdir}/train.${protein}.dose.txt
  
  echo -e "--- Training dose file generated ---\n"
  
  # Remove intermediate files
  rm ${outdir}/header.txt ${outdir}/chr*.dose.txt ${outdir}/chr*.gt.txt

  
  
  #------ Subjects and Protein level file ------#

  cut ${peddir}/WHI_olink_train.ped -f1 | sed '1d' > ${outdir}/train.subjects.txt
  
  cat ${peddir}/WHI_olink_train.ped | tr -s '\t' ',' | csvcut -c ${protein} | sed '1d' > ${outdir}/level.temp.txt
  
  paste ${outdir}/train.subjects.txt ${outdir}/level.temp.txt  > ${outdir}/train.protein_level.txt
  
  rm ${outdir}/level.temp.txt
  
  echo -e "--- Training Subjects ID file & Protein level file generated ---\n"
  
fi

# Check whether EN Training step is already done. If not, do EN Training.

if [ -f "${outdir}/train.${protein}.log" ]
then
  echo -e "--- Training step is already done. Proceed to next step. ---\n" 

else
  
  #------ EN training via R ------#
  
  Rscript "${tooldir}/pwas_training.R" \
        -g ${outdir}/train.${protein}.dose.txt \
        -p ${outdir}/train.protein_level.txt \
        -o ${outdir}/train.${protein} \
        -r 123
  
  echo -e "--- Training ${protein} by all cis + trans pQTLs EN and linear regression via R ---\n"
  
fi


echo -e "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"
echo -e "                                     Testing Step                                           \n" 
echo -e "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"

# Check whether EN model is well trained (more than 2 non-zero coefficients). If not, stop.

if grep -q "Less than 2 non-zero coefficients so no model reported" "/proj/yunligrp/users/chanhwa/pwas/whi_pwas_freeze9/result/${protein}/cisonly/train.${protein}.log"
then
  echo -e "--- EN model training failed. Stop the procedure. ---\n"
else

  # Check whether Testing data is already generated. If not, generate it.
  
  if [ -f "${outdir}/test.${protein}.dose.txt" ]
  then
    echo -e "--- Testing data is already generated. Proceed to next step. ---\n" 
  
  else
  
    #------ Testing dose file ------#
      
    start=$(( $begin - 1000000 ))
    if [[ $start -lt 1 ]]; then  
      start=1
    fi
    finish=$(( $end + 1000000 ))  
    
    if(( ${chr}==23 )); then chr=X; fi
      
    echo -e "--- Testing dose of chromosome ${chr} file generated ---\n"
    
    # Subset only cis variants from genotype data file
    bcftools view \
      -t ${chr}:$start-$finish \
      ${genodir}/test/WHI.olink.chr${chr}.dose.nochr.test.vcf.gz \
      -Oz -o ${outdir}/test.filtered.vcf.gz 
  
    # Generate header
    echo "CHROM POS REF ALT ID" $(bcftools query -l "${outdir}/test.filtered.vcf.gz") | sed 's/ /\t/g' > ${outdir}/header.txt
  
    # Extract Hard Genotype Calling (GT field) of cis variants from genotype data file
    bcftools query \
      -f '%CHROM\t%POS\t%REF\t%ALT{0}\t%ID[\t%GT]\n' \
      "${outdir}/test.filtered.vcf.gz" > ${outdir}/chr${chr}.gt.txt
    
    # Transform GT field (0/0, 0/1, 1/1) into dose (0, 1, 2)
    sed 's/0\/0/0/g' ${outdir}/chr${chr}.gt.txt | sed 's/0\/1/1/g' | sed 's/1\/1/2/g' > ${outdir}/chr${chr}.dose.txt
    
    echo -e "--- Testing dose of chromosome ${chr} file generated ---\n"
  
    # Paste header and testing dose file over chromosomes
    cat ${outdir}/header.txt ${outdir}/chr*.dose.txt > ${outdir}/test.${protein}.dose.txt
    
    echo -e "--- Testing dose file generated ---\n"
    
    # Remove intermediate files
    rm ${outdir}/header.txt ${outdir}/chr*.dose.txt ${outdir}/chr*.gt.txt
  
    
    
    #------ Subjects and Protein level file ------#
  
    cut ${peddir}/WHI_olink_test.ped -f1 | sed '1d' > ${outdir}/test.subjects.txt
    
    cat ${peddir}/WHI_olink_test.ped | tr -s '\t' ',' | csvcut -c ${protein} | sed '1d' > ${outdir}/level.temp.txt
    
    paste ${outdir}/test.subjects.txt ${outdir}/level.temp.txt  > ${outdir}/test.protein_level.txt
    
    rm ${outdir}/level.temp.txt
    
    echo -e "--- Testing Subjects ID file & Protein level file generated ---\n"
  
  fi
  
  # Check whether Testing step is already done. If not, do Testing step.
  
  if [ -f "${outdir}/test.${protein}.log" ]
  then
    echo -e "--- Testing step is already done. Proceed to next step. ---\n" 
  else
    #------ Trained model testing via R ------#    
    Rscript "${tooldir}/pwas_testing.R" \
      -g ${outdir}/test.${protein}.dose.txt \
      -p ${outdir}/test.protein_level.txt \
      -b ${outdir}/train.${protein}.betas.EN.txt \
      -l ${outdir}/train.${protein}.log \
      -o ${outdir}/test.${protein}
  
    echo -e "--- Testing trained ${protein} cis variants only EN model via R ---\n"
    
  fi

fi

echo -e "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"
echo -e "                                  All Step finished                                         \n" 
echo -e "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"

endtime=`date +%s`

runtime=$((endtime-starttime))
echo -e "Running time: $runtime seconds\n"
