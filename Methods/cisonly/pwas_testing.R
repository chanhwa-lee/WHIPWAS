#*******************************************************************************
#**********           PWAS Model Testing via Elastic Net              **********
#**********           Linear Model or Ridge Regression                **********
#**********           Used as secondary model following               **********
#**********           Variable selection via Elastic Net              **********
#**********	                                                          **********	
#**********           Written by:				                              **********
#**********           Annie Shan     - yshan@live.unc.edu             **********
#**********           Jonathan Rosen - jdrosen@live.unc.edu           **********
#**********           Munan Xie	- munan@med.unc.edu                   **********
#**********           Chanhwa Lee	- chanhwa@email.unc.edu             **********
#**********           Version: 0.7                                    **********
#**********           Nov 2, 2021                                     **********
#*******************************************************************************


# Load required libraries

suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(glmnet))
suppressPackageStartupMessages(require(data.table))
suppressPackageStartupMessages(require(dplyr))

# Set list of options for input files

option_list = list(
  make_option(c("-g", "--genotype"), action = "store", default = NA, type = "character",
              help = paste0("File name of testing genotype matrix for n subjects and p SNPs\n\t\t",
                            "Same format as the dosage file for the training step\n\t\t",
                            "A p by n+5 matrix: First five columns are\n\t\t",
                            "CHROM POS REF ALT ID(of variant)\n\t\t",
                            "and from sixth to last columns are SubjectID\n\t\t",
                            "and rows correspond to variant")),
  
  make_option(c("-p", "--protein"), action = "store", default = NA, type = "character",
              help = paste0("File name of protein level data\n\t\t",
                            "First field should be Subject ID\n\t\t",
                            "Second field should be protein level\n\t\t",
                            "adjusted for desired covariates and properly normalized")),
  
  make_option(c("-b", "--beta"), action = "store", default = NA, type = "character",
              help = paste0("File name of beta file\n\t\t",
                            "This file is the output from the training step\n\t\t",
                            "e.g. test.betas.EN.txt or test.betas.lm.txt")),
  
  make_option(c("-l", "--log"), action = "store", default = NA, type = "character",
              help = paste0("File name of log file\n\t\t",
                            "This is the output log file from pwas_training.R")),
  
  make_option(c("-o", "--output"), action = "store", default = NA, type = "character",
              help = paste0("Prefix for output file(s)\n\t\t",
                            "Suffixes include:\n\t\t",
                            ".prtex - predicted protein level\n\t\t"))
)

opt = parse_args(OptionParser(option_list = option_list))

cat(paste("*****   ",opt$o,"   *****\n") , file = paste0(opt$o, ".log"), append = T)


# Check input options

continue = TRUE
if (is.na(opt$g)) {
  cat("You must specify a genotype file\n", file = paste0(opt$o, ".log"), append = T)
  continue = FALSE
}

if (is.na(opt$p)) {
  cat("You must specify a protein file\n", file = paste0(opt$o, ".log"), append = T)
  continue = FALSE
}

if (is.na(opt$b)) {
  cat("You must specify a beta file\n", file = paste0(opt$o, ".log"), append = T)
  continue = FALSE
}

if (is.na(opt$l)) {
  cat("You must specify an output log file from the training step\n", file = paste0(opt$o, ".log"), append = T)
  continue = FALSE
}

if (is.na(opt$o)) {
  cat("You must specify a prefix for output files\n", file = paste0(opt$o, ".log"), append = T)
  continue = FALSE
}

if (! continue) {
  cat("Please correct input arguments and submit job again\n", file = paste0(opt$o, ".log"), append = T)
  quit("no") }



# #================= Example opt ====================#
# opt = data.frame(
#   g = "/proj/yunligrp/users/chanhwa/pwas/whi_pwas_freeze9/result/cardio2_ACE2/proposed/test.cardio2_ACE2.dose.txt",
#   p = "/proj/yunligrp/users/chanhwa/pwas/whi_pwas_freeze9/result/cardio2_ACE2/proposed/test.protein_level.txt",
#   b = "/proj/yunligrp/users/chanhwa/pwas/whi_pwas_freeze9/result/cardio2_ACE2/proposed/train.cardio2_ACE2.betas.EN.txt",
#   l = "/proj/yunligrp/users/chanhwa/pwas/whi_pwas_freeze9/result/cardio2_ACE2/proposed/train.cardio2_ACE2.log",
#   o = "/proj/yunligrp/users/chanhwa/pwas/whi_pwas_freeze9/result/cardio2_ACE2/proposed/test.cardio2_ACE2"
# )
# #==================================================#


## Load genotype, expression, and SNP list data

# Hard genotype values data

X = fread(opt$g, data.table=F)

# Take the hard genotype values only and make n (number of samples) by p (number of variants) matrix
X.dose = as.data.frame(t(as.matrix(X[,-(1:5)])))

# Transform characters into numbers
for(j in 1:ncol(X.dose)){ X.dose[,j] = as.numeric(X.dose[,j]) }

# "./." in the hard genotype is transformed into NA, so remove columns including NA
na.col.idx = apply(X = X.dose, FUN = function(x) any(is.na(x)), MARGIN = 2)

X.test = X.dose[,!na.col.idx]

# Let colnames of X.train as CHR:POS:REF:ALT
newID = X %>% mutate(newID = paste(CHROM, POS, REF, ALT, sep = ":")) %>% select(newID)
colnames(X.test) = newID[!na.col.idx,]

# Protein levels data
Y = read.table(opt$p, row.names = 1, header = F, stringsAsFactors = F)
y.test = Y[rownames(X.test), ]

## Check whether sample IDs perfectly overlaps
if (! (all(rownames(Y) %in% rownames(X.test)) & all(rownames(X.test) %in% rownames(Y)))) {
  cat("Subjects IDs in Genotype data and Protein data should match\n", file = paste0(opt$o, ".log"), append = T)
  cat("Please check IDs before attempting again\n\n", file = paste0(opt$o, ".log"), append = T)
  quit("no")
}


beta0     = read.table(opt$b, header = T, stringsAsFactors = F, colClasses = c("character","numeric","character","character","numeric"))
beta      = beta0[beta0[,1] != "Intercept",]
info0     = scan(opt$l, what = "char", sep = "\n", quiet = T)
info      = list()
#info$pval = as.numeric(gsub("P-value for linear model is ", "", info0[grep("P-value", info0)]))
info$rsq  = as.numeric(gsub("Elastic Net model correlation squared is ", "", info0[grep("Elastic Net model correlation", info0)]))


# Check training result
if (info$rsq < 0.01){
  cat("This prediction model is not very good in terms of p-value or adjusted r-squared, so we won't use it for predicting prtex.\n", file = paste0(opt$o, ".log"), append = T)
  quit("no")
}

# Get counts for overlapping SNPs
rownames(beta) = paste(beta[,1], beta[,2], beta[,3], beta[,4], sep = ":")
n.snps.dosage = length(X.test)
n.snp.beta = dim(beta)[1]
all.snps.overlap = all(rownames(beta) %in% colnames(X.test))

if (!all.snps.overlap) {
  cat("Some variants in the prediction model are not found in the input genotype file,\n", file = paste0(opt$o, ".log"), append = T)
  cat("Restrict variants only in the input genotype file\n", file = paste0(opt$o, ".log"), append = T)
  beta = beta[rownames(beta) %in% colnames(X.test),]
}

cat(paste0("Number of variants used in prediction is ", nrow(beta), "\n"), file = paste0(opt$o, ".log"), append = T)

# Predict protein level

X.test.trim = as.matrix(X.test[ , match(rownames(beta), colnames(X.test))])

prtlvl = data.frame(predicted = X.test.trim%*%beta[,5], observed = y.test)

write.table(prtlvl, file = paste0(opt$o, ".prtex.PearsonR2=",cor(prtlvl$predicted, prtlvl$observed)**2,".txt"), row.names = T, col.names = T, quote = F)
cat("Pearson R2 of predicted and observed protein level is", cor(prtlvl$predicted, prtlvl$observed)**2,"\n", file = paste0(opt$o, ".log"), append = T)
