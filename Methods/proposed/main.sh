#!/bin/bash

#----- Variables for proposed method -----#

starttime=`date +%s`

protein=$1
pval_thresh=$2

genodir="/proj/yunligrp/users/chanhwa/pwas/whi_pwas_freeze9/data/genotypedata"

peddir="/proj/yunligrp/users/chanhwa/pwas/whi_pwas_freeze9/data/proteindata"

tooldir="/proj/yunligrp/users/chanhwa/pwas/whi_pwas_freeze9/method/proposed"

outdir="/proj/yunligrp/users/chanhwa/pwas/whi_pwas_freeze9/result/${protein}/proposed_pval_${pval_thresh}"

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

# Check whether kinship matrix is already generated. If not, generate it.

if [ -f ${tooldir}/WHI.olink.kinf ]
then
  echo -e "--- Kinship matrix is already generated. Proceed to next step. ---\n" 
else
  epacts make-kin \
    --vcf ${genodir}/all/WHI.olink.chr1.dose.nochr.all.vcf.gz \
    --sepchr \
    --out ${tooldir}/WHI.olink.kinf  \
    --run 5
fi

# Check whether GWAS step is already done. If not, do it.

cp -R "/proj/yunligrp/users/chanhwa/pwas/whi_pwas_freeze9/result/${protein}/proposed/gwas" ${outdir}/gwas

if [ -f ${outdir}/gwas/${protein}_gwas.epacts.gz ]
then
  echo -e "--- GWAS step is already done. Proceed to next step. ---\n" 
else
  mkdir -p ${outdir}/gwas

  epacts single \
    --vcf ${genodir}/gwas/WHI.olink.chr1.dose.nochr.gwas.vcf.gz \
    --sepchr \
    --ped ${peddir}/WHI_olink_gwas.ped \
    --pheno ${protein} \
    --cov plate_1 --cov plate_2 --cov plate_3 --cov plate_4 --cov plate_5 \
    --cov plate_6 --cov plate_7 --cov plate_8 --cov plate_9 --cov plate_10 \
    --cov plate_11 --cov plate_12 --cov plate_13 --cov plate_14 --cov plate_15 \
    --min-maf 0.01 \
    --kin ${tooldir}/WHI.olink.kinf \
    --test q.emmax \
    --out ${outdir}/gwas/${protein}_gwas \
    --run 5
    
fi


# Check whether Training data is already generated. If not, generate it.

if [ -f "${outdir}/train.${protein}.dose.txt" ]
then
  echo -e "--- Training data is already generated. Proceed to next step. ---\n" 

else


  #------ GWAS significant pQTLs ------#
  
  zcat ${outdir}/gwas/${protein}_gwas.epacts.gz | awk -v pval=${pval_thresh} '$11 < pval {print $4}' | cut -d '_' -f3 > ${outdir}/pQTLs.txt
  
  echo -e "--- GWAS significant pQTLs list generated. ---\n"


  #------ Training dose file ------#
  
  for chr in {1..23};do
  
    if(( ${chr}==23 )); then chr=X; fi
  
    # Extract Hard Genotype Calling (GT field) of cis variants from genotype data file
    bcftools query \
      -i ID=@${outdir}/pQTLs.txt \
      -f '%CHROM\t%POS\t%REF\t%ALT{0}\t%ID[\t%GT]\n' \
      ${genodir}/pwas/WHI.olink.chr${chr}.dose.nochr.pwas.vcf.gz > ${outdir}/chr${chr}.gt.txt
    
    # Transform GT field (0/0, 0/1, 1/1) into dose (0, 1, 2)
    sed 's/0\/0/0/g' ${outdir}/chr${chr}.gt.txt | sed 's/0\/1/1/g' | sed 's/1\/1/2/g' > ${outdir}/chr${chr}.dose.txt
    
    echo -e "--- Training dose of chromosome ${chr} file generated ---\n"
  
  done
  
  # Generate header
  echo "CHROM POS REF ALT ID" $(bcftools query -l ${genodir}/pwas/WHI.olink.chr${chr}.dose.nochr.pwas.vcf.gz) | sed 's/ /\t/g' > ${outdir}/header.txt
  
  # Paste header and training dose file over chromosomes
  cat ${outdir}/header.txt ${outdir}/chr*.dose.txt > ${outdir}/train.${protein}.dose.txt
  
  echo -e "--- Training dose file generated ---\n"
  
  # Remove intermediate files
  rm ${outdir}/header.txt ${outdir}/chr*.dose.txt ${outdir}/chr*.gt.txt

  
  
  #------ Subjects and Protein level file ------#

  cut ${peddir}/WHI_olink_pwas.ped -f1 | sed '1d' > ${outdir}/train.subjects.txt
  
  cat ${peddir}/WHI_olink_pwas.ped | tr -s '\t' ',' | csvcut -c ${protein} | sed '1d' > ${outdir}/level.temp.txt
  
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

if grep -q "Less than 2 non-zero coefficients so no model reported" "${outdir}/train.${protein}.log"
then
  echo -e "--- EN model training failed. Stop the procedure. ---\n"
else

  # Check whether Testing data is already generated. If not, generate it.

  if [ -f "${outdir}/test.${protein}.dose.txt" ]
  then
    echo -e "--- Testing data is already generated. Proceed to next step. ---\n" 
  else
  
    #------ Testing dose file ------#
    
    for chr in {1..23};do
    
      if(( ${chr}==23 )); then chr=X; fi
    
      # Extract Hard Genotype Calling (GT field) of cis variants from genotype data file
      bcftools query \
        -i ID=@${outdir}/pQTLs.txt \
        -f '%CHROM\t%POS\t%REF\t%ALT{0}\t%ID[\t%GT]\n' \
        ${genodir}/test/WHI.olink.chr${chr}.dose.nochr.test.vcf.gz > ${outdir}/chr${chr}.gt.txt
      
      # Transform GT field (0/0, 0/1, 1/1) into dose (0, 1, 2)
      sed 's/0\/0/0/g' ${outdir}/chr${chr}.gt.txt | sed 's/0\/1/1/g' | sed 's/1\/1/2/g' > ${outdir}/chr${chr}.dose.txt
      
      echo -e "--- Testing dose of chromosome ${chr} file generated ---\n"
    
    done
    
    # Generate header
    echo "CHROM POS REF ALT ID" $(bcftools query -l ${genodir}/test/WHI.olink.chr${chr}.dose.nochr.test.vcf.gz) | sed 's/ /\t/g' > ${outdir}/header.txt
    
    # Paste header and training dose file over chromosomes
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
