#*******************************************************************************
#**********           PWAS Model Training via Elastic Net             **********
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
#**********           Dec 2, 2020                                     **********
#*******************************************************************************


# Load required libraries

suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(glmnet))
suppressPackageStartupMessages(require(data.table))
suppressPackageStartupMessages(require(R.utils))

# Set list of options for input files

option_list = list(
  make_option(c("-d", "--dosage"), action = "store", default = NA, type = "character",
              help = paste0("File name of genotype matrix for n subjects and p SNPs\n\t\t",
                            "Same format as the dosage file for the training step\n\t\t",
                            "Format for SNP ID should be chr#:x, where # is chromosome number\n\t\t",
                            "And x is the chromosome position under BUILD 38\n\t\t",
                            "The subsequent n lines should contain Subject ID in field 1\n\t\t",
                            "and dosages in fields 2:(p+1)")),
  make_option(c("-e", "--expression"), action = "store", default = NA, type = "character",
              help = paste0("File name of expression(protein) data\n\t\t",
                            "First field should be Subject ID\n\t\t",
                            "Second field should be expression value\n\t\t",
                            "adjusted for desired covariates and properly normalized")),
  make_option(c("-s", "--snps"), action = "store", default = NA, type = "character",
              help = paste0("File name of list of SNPs\n\t\t",
                            "This file should be four fields with one SNP ID per line\n\t\t",
                            "SNP ID format is chr#, where # is chromosome number\n\t\t",
                            "Next field is the chromosome position under BUILD 38\n\t\t",
                            "The next two fields should be dosage allele and other allele\n\t\t",
                            "In that order")),
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
                            ".prtex - predicted protein expression level\n\t\t",
                            ".prtex2 - predicted protein expression level scaled to [0,2]"))
)

opt = parse_args(OptionParser(option_list = option_list))


# Check input options

continue = TRUE
if (is.na(opt$d)) {
  cat("You must specify a dosage file\n")
  continue = FALSE
}

if (is.na(opt$s)) {
  cat("You must specify a snplist file\n")
  continue = FALSE
}

if (is.na(opt$b)) {
  cat("You must specify a beta file\n")
  continue = FALSE
}

if (is.na(opt$l)) {
  cat("You must specify an output log file from the training step\n")
  continue = FALSE
}

if (is.na(opt$o)) {
  cat("You must specify a prefix for output files\n")
  continue = FALSE
}

if (! continue) {
  cat("Please correct input arguments and submit job again\n")
  quit("no") }

# ################## TEST #################################
# opt <- data.frame(d="WHI_olink.neurol_MANF.dosage.gz",
#                   e="WHI_olink.neurol_MANF.expression.txt",
#                   s="WHI_olink.neurol_MANF.snps",
#                   b="WHI_olink.neurol_MANF.betas.EN.txt",
#                   l="WHI_olink.neurol_MANF.log",
#                   o="WHI_olink.neurol_MANF", stringsAsFactors = F)
# 
# opt
# #########################################################


# Load genotype, expression, and SNP list data

input.file = file(opt$d)

if (summary(input.file)$class == "gzfile") {
        dosage0 = fread(opt$d, data.table = F, header = F)
} else {
        dosage0 = fread(input.file, data.table = F, header = F)
        cat("Dosage file is not zipped!\n")
        cat("Save storage space... compress your data!\n")
}

close(input.file)

row.names(dosage0) = dosage0[,1]
dosage0 = dosage0[, -1]

snp.list = read.table(opt$s, header = F, stringsAsFactors = F, colClasses = "character")

beta0     = read.table(opt$b, header = T, stringsAsFactors = F, colClasses = c("character","numeric","character","character","numeric"))
beta      = beta0[beta0[,1] != "Intercept",]
info0     = scan(opt$l, what = "char", sep = "\n", quiet = T)
info      = list()
#info$pval = as.numeric(gsub("P-value for linear model is ", "", info0[grep("P-value", info0)]))
info$rsq  = as.numeric(gsub("Elastic Net model correlation squared is ", "", info0[grep("Elastic Net model correlation", info0)]))

# Check dosage file

if (any(dosage0 > 2 || dosage0 < 0 || is.na(dosage0))) {
  cat("Dosage values must be between 0 and 2 inclusive and can not contain missing values\n")
  cat("Please check values before attempting again\n")
  quit("no")
}


if (anyDuplicated(rownames(dosage0))) {
  cat("Found duplicated subject IDs in dosage file\n")
  cat("Please check file before attempting again\n")
  quit("no")
}

# Check SNP list file

if (dim(snp.list)[2] != 4) {
        cat("SNP list file must have four fields\n")
        cat("Format is chr# pos allele1 allele2\n")
        cat("Check file before attempting again\n\n")
        quit("no")
}

## Remove duplicate SNPs (due to all cis + trans pQTLs)

snp.id = paste(snp.list[,1], snp.list[,2], snp.list[,3], snp.list[,4], sep = ":")
uniq.snp.id = unique(snp.id)
uniq.index = match(uniq.snp.id, snp.id)
snp.list = snp.list[uniq.index,]
dosage0 = dosage0[,uniq.index]

rownames(snp.list) = paste(snp.list[,1], snp.list[,2], snp.list[,3], snp.list[,4], sep = ":")
colnames(dosage0) = rownames(snp.list)

# Check training result
if (info$rsq < 0.01){
  cat("This prediction model is not very good in terms of p-value or adjusted r-squared, so we won't use it for predicting prtex.\n")
  quit("no")
}

# Get counts for overlapping SNPs
rownames(beta) = paste(beta[,1], beta[,2], beta[,3], beta[,4], sep = ":")
n.snps.dosage = length(dosage0)
n.snp.beta = dim(beta)[1]
all.snps.overlap = all(rownames(beta) %in% colnames(dosage0))

if (!all.snps.overlap) {
  cat("Some SNPs in the prediction model are not found in the input dosage file, exit.\n")
  quit("no")
}


observed.prt = read.table(opt$e, header = F, stringsAsFactors = F, colClasses = c("character","numeric"))
row.names(observed.prt) = observed.prt[,1]

dosage = as.matrix(dosage0[match(observed.prt[,1], rownames(dosage0)),match(rownames(beta), colnames(dosage0))])

# predict prtex (protein expression level)
prtex = data.frame(predicted = dosage%*%beta[,5], observed = observed.prt[,2])

write.table(prtex, file = paste0(opt$o, ".prtex.PearsonR2=",cor(prtex$predicted, prtex$observed)**2,".txt"), row.names = T, col.names = T, quote = F)
cat("Pearson R2 of predicted and observed protein level is", cor(prtex$predicted, prtex$observed)**2,"\n")



#######################################################
# # scale predicted prtex to [0,2]
# prtex0 = prtex - min(prtex)
# prtex2 = 2*prtex0/max(prtex0)
# write.table(prtex2, file = paste0(opt$o, ".prtex2"), row.names = T, col.names = F, quote = F)