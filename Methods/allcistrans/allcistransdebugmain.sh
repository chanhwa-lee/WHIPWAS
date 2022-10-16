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