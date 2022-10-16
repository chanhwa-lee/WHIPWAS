#!/bin/bash

starttime=`date +%s`

protein=$1
chr=$2
begin=$3
end=$4

proteindir="/proj/yunligrp/users/chanhwa/pwas/whi_pwas_neat/result/${protein}"

mkdir -p ${proteindir}

outdir="${proteindir}/allcistrans"

mkdir -p ${outdir}

cd ${outdir}


echo
echo "----------------------------------------------------------------------------------------------"
echo "----------------------------------------------------------------------------------------------"
echo ""
echo "          Protein ${protein} all Cis and trans pQTLs prediction EN model training and testing                " 
echo ""
echo "                              `date`                                                          "
echo ""
echo "----------------------------------------------------------------------------------------------"
echo "----------------------------------------------------------------------------------------------"
echo

#------ Load Modules ------#

module add python
module add samtools
module add r/4.0.3


echo
echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
echo ""
echo "                                       Training Step                                          " 
echo ""
echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
echo


#------ GWAS significant pQTLs ------#

pval_thresh=0.0001

zcat ${proteindir}/gwas/${protein}_gwas.epacts.gz | awk -v pval=${pval_thresh} '$11 < pval {print $4}' | cut -d '_' -f3 > ${outdir}/pQTLs.txt

echo "--- GWAS significant pQTLs list generated ---"
echo

#------ Subjects and Protein level file ------#

ped_file="/proj/yunligrp/users/chanhwa/pwas/whi_pwas_neat/1_dataset/proteindata/WHI_olink_pwas.ped"

cut ${ped_file} -f1 | sed '1d' > ${outdir}/train.subjects.txt

cat ${ped_file} | tr -s '\t' ',' | csvcut -c ${protein} | sed '1d' > ${outdir}/level.temp.txt

paste ${outdir}/train.subjects.txt ${outdir}/level.temp.txt  > ${outdir}/train.protein_level.txt

rm ${outdir}/level.temp.txt

echo "--- Training Subjects ID file & Protein level file generated ---"
echo

#------ Dosage Conversion for EN model training input ------#

##----------------- All cis information -------------------##

genodir="/proj/yunligrp/users/chanhwa/pwas/whi_pwas_neat/1_dataset/genotype"

start=$(( $begin - 1000000 ))
if [[ $start -lt 1 ]]; then  
  start=1
fi
finish=$(( $end + 1000000 ))  

if(( ${chr}==23 )); then chr=X; fi
  
echo "*** Training Dosage Converting for chr ${chr} (cis variants) ***"
              
bcftools view \
  -t chr${chr}:$start-$finish \
  ${genodir}/proposed_train/WHI.olink.chr${chr}.dose.proposed_train.vcf.gz \
  -Oz -o ${outdir}/train.filtered.vcf.gz 

/proj/yunligrp/users/chanhwa/pwas/whi_pwas_neat/allcistrans/DosageConvertor \
  --vcfDose ${outdir}/train.filtered.vcf.gz \
  --prefix ${outdir}/train.${protein} \
  --type mach \
  --format 1

rm train.filtered.vcf.gz
zcat train.${protein}.mach.dose.gz | cut -f 1,3- | cut -d '>' -f2- > train.${protein}.dosage
rm train.${protein}.mach.dose.gz
cat train.${protein}.mach.info | sed '1d' | cut -f1 | sed 's/:/\t/g' > train.${protein}.snps

echo
echo "--- Making <train.${protein}.dosage> and <train.${protein}.snps> files for all cis variants for training ---"
echo

##----------------- trans pQTLs information -------------------##

for chrom in {1..23}; do

  if(( ${chrom}==23 )); then chrom=X; fi
  
  echo "*** Training Dosage Converting for chr ${chrom} (trans pQTLs) ***"
  
  bcftools view -i ID=@${outdir}/pQTLs.txt \
      ${genodir}/proposed_train/WHI.olink.chr${chrom}.dose.proposed_train.vcf.gz \
      -Oz -o ${outdir}/train.chr${chrom}.filtered.vcf.gz
      
  /proj/yunligrp/users/chanhwa/pwas/whi_pwas_neat/allcistrans/DosageConvertor \
      --vcfDose ${outdir}/train.chr${chrom}.filtered.vcf.gz \
      --prefix ${outdir}/train.chr${chrom}.${protein} \
      --type mach \
      --format 1
done

###--------- Merging converted dosage files (cis & trans) ---------###

for chrom in {1..23}; do

  if(( ${chrom}==23 )); then chrom=X; fi
  
  zcat train.chr${chrom}.${protein}.mach.dose.gz | cut -f3- > temp
  paste train.${protein}.dosage temp > new_train.${protein}.dosage
  mv new_train.${protein}.dosage train.${protein}.dosage
  
  cat train.chr${chrom}.${protein}.mach.info | sed '1d' | cut -f1 | sed 's/:/\t/g' >> train.${protein}.snps
  
done

cat train.${protein}.dosage | bgzip > train.${protein}.dosage.gz

echo
echo "--- Making <train.${protein}.dosage.gz> and <train.${protein}.snps> files for training ---"
echo


#------ EN training via R ------#

Rscript /proj/yunligrp/users/chanhwa/pwas/whi_pwas_neat/allcistrans/allcistrans_training.R \
        -d ${outdir}/train.${protein}.dosage.gz \
        -e ${outdir}/train.protein_level.txt \
        -s ${outdir}/train.${protein}.snps \
        -o ${outdir}/train.${protein} \
        -r 123

echo
echo "--- Training ${protein} by all cis + trans pQTLs EN and linear regression via R ---"
echo


#------ Clean Directory ------#

mkdir -p mediate_files
mv ${outdir}new_train.${protein}.dosage ${outdir}/temp ${outdir}/train.chr*.${protein}.mach.* ${outdir}/*_dosage_convert.log ${outdir}/train.chr*.filtered.vcf.gz mediate_files


echo
echo
echo
echo
echo
echo
echo
echo
echo
echo
echo
echo
echo
echo
echo
echo
echo
echo
echo
echo
echo
echo
echo
echo


echo
echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
echo ""
echo "                                       Testing Step                                           " 
echo ""
echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
echo

#------ GWAS significant pQTLs ------#

if [[ ! -f ${outdir}/pQTLs.txt ]]
then
    echo "pQTL list file does not exist. Check it."
fi

echo

#------ Testing Subjects and Protein level file ------#

ped_file="/proj/yunligrp/users/chanhwa/pwas/whi_pwas_neat/1_dataset/proteindata/WHI_olink_test.ped"

cut ${ped_file} -f1 | sed '1d' > ${outdir}/test.subjects.txt

cat ${ped_file} | tr -s '\t' ',' | csvcut -c ${protein} | sed '1d' > ${outdir}/level.temp.txt

paste ${outdir}/test.subjects.txt ${outdir}/level.temp.txt  > ${outdir}/test.protein_level.txt

rm ${outdir}/level.temp.txt

echo "--- Testing Subjects ID file & Protein level file generated ---"
echo

#------ Dosage Conversion for EN model testing input ------#

##----------------- All cis information -------------------##
  
echo "*** Testing Dosage Converting for chr ${chr} (cis variants) ***"
              
bcftools view \
  -t chr${chr}:$start-$finish \
  ${genodir}/general_test/WHI.olink.chr${chr}.dose.general_test.vcf.gz \
  -Oz -o ${outdir}/test.filtered.vcf.gz 

/proj/yunligrp/users/chanhwa/pwas/whi_pwas_neat/allcistrans/DosageConvertor \
  --vcfDose ${outdir}/test.filtered.vcf.gz \
  --prefix ${outdir}/test.${protein} \
  --type mach \
  --format 1

rm test.filtered.vcf.gz
zcat test.${protein}.mach.dose.gz | cut -f 1,3- | cut -d '>' -f2- > test.${protein}.dosage
rm test.${protein}.mach.dose.gz
cat test.${protein}.mach.info | sed '1d' | cut -f1 | sed 's/:/\t/g' > test.${protein}.snps

echo
echo "--- Making <test.${protein}.dosage.gz> and <test.${protein}.snps> files for testing ---"
echo

##----------------- trans pQTLs information -------------------##

for chrom in {1..23}; do

  if(( ${chrom}==23 )); then chrom=X; fi
  
  echo "*** Testing Dosage Converting for chr ${chrom} (trans pQTLs) ***"
  
  bcftools view -i ID=@${outdir}/pQTLs.txt \
      ${genodir}/general_test/WHI.olink.chr${chrom}.dose.general_test.vcf.gz \
      -Oz -o ${outdir}/test.chr${chrom}.filtered.vcf.gz
  
  /proj/yunligrp/users/chanhwa/pwas/whi_pwas_neat/allcistrans/DosageConvertor \
      --vcfDose ${outdir}/test.chr${chrom}.filtered.vcf.gz \
      --prefix ${outdir}/test.chr${chrom}.${protein} \
      --type mach \
      --format 1
done


###--------- Merging converted dosage files (cis & trans) ---------###

for chrom in {1..23}; do

  if(( ${chrom}==23 )); then chrom=X; fi
  
  zcat test.chr${chrom}.${protein}.mach.dose.gz | cut -f3- > temp
  paste test.${protein}.dosage temp > new_test.${protein}.dosage
  mv new_test.${protein}.dosage test.${protein}.dosage
  
  cat test.chr${chrom}.${protein}.mach.info | sed '1d' | cut -f1 | sed 's/:/\t/g' >> test.${protein}.snps
  
done

cat test.${protein}.dosage | bgzip > test.${protein}.dosage.gz

echo
echo "--- Making <test.${protein}.dosage.gz> and <test.${protein}.snps> files for testing ---"
echo


#------ Trained model testing via R ------#

Rscript /proj/yunligrp/users/chanhwa/pwas/whi_pwas_neat/allcistrans/allcistrans_prediction.R \
  -d ${outdir}/test.${protein}.dosage.gz \
  -e ${outdir}/test.protein_level.txt \
  -s ${outdir}/test.${protein}.snps \
  -b ${outdir}/train.${protein}.betas.EN.txt \
  -l ${outdir}/train.${protein}.log \
  -o ${outdir}/test.${protein}

echo
echo "--- Testing trained ${protein} all cis + trans pQTLs EN model via R ---"
echo

#------ Clean Directory ------#

mkdir -p mediate_files
mv ${outdir}/new_test.${protein}.dosage ${outdir}/temp ${outdir}/test.chr*.${protein}.mach.* ${outdir}/test.chr*.filtered.vcf.gz mediate_files


echo
echo "**********************************************************************************************"
echo ""
echo "                                      All step finished                                       " 
echo ""
echo "**********************************************************************************************"
echo

endtime=`date +%s`

runtime=$((endtime-starttime))
echo "Running time: $runtime seconds"