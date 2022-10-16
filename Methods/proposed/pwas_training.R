#*******************************************************************************
#**********           PWAS Model Training via Elastic Net             **********
#**********           Linear Model or Ridge Regression                **********
#**********           Used as secondary model following               **********
#**********           Variable selection via Elastic Net              **********
#**********                                                           **********	
#**********           Written by:				                              **********
#**********           Jonathan Rosen - jdrosen@live.unc.edu           **********
#**********           Annie Shan - yshan@live.unc.edu                 **********
#**********           revised by Munan on March19, 2019               **********
#**********           revised by Chanhwa on Oct26, 2021               **********
#**********           Version: 0.8                                    **********
#**********           Oct 26, 2021                                    **********
#*******************************************************************************


## Load required libraries

suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(glmnet))
suppressPackageStartupMessages(require(data.table))
suppressPackageStartupMessages(require(dplyr))
if(!require(glmnetUtils))  install.packages("glmnetUtils")
suppressPackageStartupMessages(require(glmnetUtils))

## Start timer for code

time = proc.time()


## Set list of options for input files

option_list = list(
  make_option(c("-g", "--genotype"), action = "store", default = NA, type = "character",
              help = paste0("File name of training genotype matrix for n subjects and p SNPs\n\t\t",
                            "A p by n+5 matrix: First five columns are\n\t\t",
                            "CHROM POS REF ALT ID(of variant)\n\t\t",
                            "and from sixth to last columns are SubjectID\n\t\t",
                            "and rows correspond to variant")),
  make_option(c("-p", "--protein"), action = "store", default = NA, type = "character",
              help = paste0("File name of training protein level data\n\t\t",
                            "First field should be Subject ID\n\t\t",
                            "Second field should be protein level value\n\t\t",
                            "adjusted for desired covariates and properly normalized")),
  make_option(c("-r", "--randomseed"), action = "store", default = NA, type = "integer",
              help = paste0("A number that used for set.seed")),
  make_option(c("-o", "--output"), action = "store", default = NA, type = "character",
              help = paste0("Prefix for output file(s)\n\t\t",
                            "Suffixes include:\n\t\t",
                            ".betas.EN.txt - elastic net model with alpha = 0.5\n\t\t",
                            ".betas.lm.txt - linear model fit using non-zero coefficients from EN\n\t\t",
                            ".betas.ridge.txt - ridge regression in case lm not fit due to collinearity"))
)

opt = parse_args(OptionParser(option_list = option_list))

#================= Example opt ====================#
opt = data.frame(
  g = "/proj/yunligrp/users/chanhwa/pwas/whi_pwas_freeze9/result/cardio2_ACE2/proposed_pval_0.0005/train.cardio2_ACE2.dose.txt",
  p = "/proj/yunligrp/users/chanhwa/pwas/whi_pwas_freeze9/result/cardio2_ACE2/proposed_pval_0.0005/train.protein_level.txt",
  r = 123,
  o = "/proj/yunligrp/users/chanhwa/pwas/whi_pwas_freeze9/result/cardio2_ACE2/proposed_pval_0.0005/example.train.protein"
)
#==================================================#

## Write output directory to stdout

cat(paste("*****   ",opt$o,"   *****\n"))


## Load genotype and protein level data

# Hard genotype values data

X = fread(opt$g, data.table=F)

# Take the hard genotype values only and make n (number of samples) by p (number of variants) matrix
X.dose = as.data.frame(t(as.matrix(X[,-(1:5)])))

# Transform characters into numbers
for(j in 1:ncol(X.dose)){ X.dose[,j] = as.numeric(X.dose[,j]) }

# "./." in the hard genotype is transformed into NA, so remove columns including NA
na.col.idx = apply(X = X.dose, FUN = function(x) any(is.na(x)), MARGIN = 2)

X.train = X.dose[,!na.col.idx]

# Let colnames of X.train as CHR:POS:REF:ALT
newID = X %>% mutate(newID = paste(CHROM, POS, REF, ALT, sep = ":")) %>% select(newID)
colnames(X.train) = newID[!na.col.idx,]

# Protein levels data
Y = read.table(opt$p, row.names = 1, header = F, stringsAsFactors = F)
y.train = Y[rownames(X.train), ]

## Check whether sample IDs perfectly overlaps
if (! (all(rownames(Y) %in% rownames(X.train)) & all(rownames(X.train) %in% rownames(Y)))) {
  cat("Subjects IDs in Genotype data and Protein data should match\n")
  cat("Please check IDs before attempting again\n\n")
  quit("no")
}

## Check genotype file

if (!all(as.matrix(X.train) %in% c(0,1,2))) {
  cat("Genotype values must be 0 or 1 or 2 and can not contain missing values\n")
  cat("Please check values before attempting again\n\n")
  quit("no")
}

if (anyDuplicated(rownames(X.train))) {
  cat("Found duplicated subject IDs in dosage file\n")
  cat("Please check file before attempting again\n\n")
  quit("no")
}


## Check protein level

if (! is.numeric(y.train) || any(is.na(y.train))) {
  cat("Protein levels must be numeric and non-missing\n")
  cat("Please check values before attempting again\n\n")
  quit("no")
}

if (anyDuplicated(rownames(y.train))) {
  cat("Found duplicated subject IDs in protein level file\n")
  cat("Please check file before attempting again\n\n")
  quit("no")
}


## Rectify counts for SNPs

n.snps.geno = dim(X.train)[2]

cat(paste0("There are ", n.snps.geno, " SNPs in the dosage file\n"))

cat("************ SNP Information ************\n", file = paste0(opt$o, ".log"), append = F)
cat(paste0("There are ", n.snps.geno, " SNPs in the genotype file\n"), file = paste0(opt$o, ".log"), append = T)

## Get counts for overlapping subjects

n.subj.geno = length(rownames(X.train))
n.subj.prot = length(y.train)

cat(paste0("There are ", n.subj.geno, " subjects in the genotype file\n"))
cat(paste0("There are ", n.subj.prot, " subjects in the protein file\n"))

cat("\n\n************ Subject Information ************\n", file = paste0(opt$o, ".log"), append = T)
cat(paste0("There are ", n.subj.geno, " subjects in the genotype file\n"), file = paste0(opt$o, ".log"), append = T)
cat(paste0("There are ", n.subj.prot, " subjects in the protein file\n"), file = paste0(opt$o, ".log"), append = T)


## Test normality of protein values

normality.test = shapiro.test(y.train)
shapiro.res = prettyNum(normality.test$p.val, format = "fg", digits = 5)

if (shapiro.res < 0.05) {
  cat(paste0("Warning: Shapiro-Wilk test for normality of expression values resulted in p-value of ", shapiro.res, "\n"))
  cat("Analysis will continue but it is recommended to check expression values\n")
  cat("\n\n************ Test for normality of expression values ************\n", file = paste0(opt$o, ".log"), append = T)
  cat(paste0("Warning: Shapiro-Wilk test for normality of expression values resulted in p-value of ", shapiro.res, "\n"), file = paste0(opt$o, ".log"), append = T)
  cat("Analysis will continue but it is recommended to check expression values\n", file = paste0(opt$o, ".log"), append = T)
} else {
  cat(paste0("Shapiro-Wilk test for normality of expression values resulted in p-value of ", shapiro.res, "\n"))
  cat("\n\n************ Test for normality of expression values ************\n", file = paste0(opt$o, ".log"), append = T)
  cat(paste0("Shapiro-Wilk test for normality of expression values resulted in p-value of ", shapiro.res, "\n"), file = paste0(opt$o, ".log"), append = T)
}


## remove perfectly correlated SNPs

trim.dosage = function(X) {
  rsq = cor(X)
  rsq.rm = (abs(rsq - diag(ncol(X))) == 1)
  snps = which(apply(rsq.rm, 1, any))
  pool = snps
  keep = NULL
  dif = setdiff(pool, keep)
  
  while(length(dif) > 0){
    keep = c(keep, dif[1])
    rm = setdiff(which(rsq.rm[dif[1],]), dif[1])
    pool = setdiff(pool, rm)
    dif = setdiff(pool, keep)
  }
  
  snps.rm = snps[! snps %in% keep]
  
  if(length(snps.rm) > 0) {
    X.trim = X[, -snps.rm]
    return(X.trim)
  } else {
    return(X)
  }
}


X.train.trim = trim.dosage(X.train)

n.trimmed = dim(X.train)[2] - dim(X.train.trim)[2]

cat(paste0("There were ", n.trimmed, " SNPs removed due to perfect collinearity\n"))
cat("\n\n************ Collinearity Information ************\n", file = paste0(opt$o, ".log"), append = T)
cat(paste0("There were ", n.trimmed, " SNPs removed due to perfect collinearity\n"), file = paste0(opt$o, ".log"), append = T)


#==============================================================================#
#==============================================================================#
#==============================================================================#

# EN Goal: minimize MSE + lambda*(alpha * L1 + (1-alpha) * L2)

## 1. Fit Elastic Net 

set.seed(opt$r)
fit.tuning = cva.glmnet(as.matrix(X.train.trim), y.train, family = "gaussian", standardize = F)

# Get optimal alpha hyper parameter.
alpha <- fit.tuning$alpha
lambdaMin <- sapply(fit.tuning$modlist, `[[`, "lambda.min")
lambdaSE <- sapply(fit.tuning$modlist, `[[`, "lambda.1se")
error <- sapply(fit.tuning$modlist, function(mod) {min(mod$cvm)})
best <- which.min(error)
optimal.tuning = data.frame(alpha = alpha[best], lambdaMin = lambdaMin[best],
                            lambdaSE = lambdaSE[best], eror = error[best])

cat(paste0("Optimal alpha hyper parameter for the Elastic Net model is ", optimal.tuning$alpha, "\n"))
cat("\n\n************ Elastic Net Model Information ************\n", file = paste0(opt$o, ".log"), append = T)
cat(paste0("Optimal alpha hyper parameter for the Elastic Net model is ", optimal.tuning$alpha, "\n"), file = paste0(opt$o, ".log"), append = T)


# Model fitting with optimal alpha hyper parameter.
fit = cv.glmnet(as.matrix(X.train.trim), y.train, 
                alpha = optimal.tuning$alpha, family = "gaussian", standardize = F)

fit.df = data.frame(fit$cvm, fit$lambda, 1:length(fit$cvm))

best.lam = fit.df[fit.df[,2] == fit$lambda.min, ]
ret = as.data.frame(fit$glmnet.fit$beta[, best.lam[1,3]])

if(all(ret[,1] == 0)){
  cat("Best model coefficient is estimated by zero vector - no SNPs are predictive \n")
} else if (var(as.matrix(X.train.trim) %*% ret[,1]) == 0) {
  cat("Model is singular - no correlation reported\n")
} else {
  model.corsq = prettyNum(cor(as.matrix(X.train.trim) %*% ret[,1], y.train)^2, format = "fg", digits = 5)
  cat(paste0("Elastic Net model correlation squared is ", model.corsq, "\n"))
  cat(paste0("Elastic Net model correlation squared is ", model.corsq, "\n"), file = paste0(opt$o, ".log"), append = T)
}

beta.values = ret[ret[,1] != 0, 1]
n.betas = length(beta.values)
if (n.betas < 2) {
  cat("Less than 2 non-zero coefficients so no model reported\n\n")
  cat("Less than 2 non-zero coefficients so no model reported\n", file = paste0(opt$o, ".log"), append = T)
  quit("no")
}

beta.snps = colnames(X.train.trim)[ret[,1] != 0]
beta.snps.split = matrix(unlist(strsplit(beta.snps, split = ":")), nrow = length(beta.snps), ncol = 4, byrow = T)
betas.elnet = cbind(beta.snps.split, beta.values)
colnames(betas.elnet) = c("Chr", "Position", "Dosage_allele", "Other_allele", "Effect")
n.model = dim(X.train.trim)[2]

cat(paste0("There are ", n.betas, " non-zero coefficients out of ", n.model, " SNPs used to fit the model\n"))
cat(paste0("There are ", n.betas, " non-zero coefficients out of ", n.model, " SNPs used to fit the model\n"), file = paste0(opt$o, ".log"), append = T)

write.table(betas.elnet, file = paste0(opt$o, ".betas.EN.txt"), col.names = T, row.names = F, quote = F, sep = "\t")

#==============================================================================#
#==============================================================================#
#==============================================================================#



## 2. Use non-zero coefficient SNPs and fit linear model

data = as.data.frame(cbind(y.train, X.train.trim[, beta.snps]))
colnames(data) = c("y", colnames(X.train.trim[, beta.snps]))


lm.summary = summary(lm(y ~ ., data = data))
if (any(lm.summary$aliased)) {
  cat("Singularities in the attempted linear model\n")
  cat("Fitting ridge regression model instead\n")
  cat("\n\n************ Secondary Model Information ************\n", file = paste0(opt$o, ".log"), append = T)
  cat("Singularities in the attempted linear model\n", file = paste0(opt$o, ".log"), append = T)
  cat("Fitting ridge regression model instead\n", file = paste0(opt$o, ".log"), append = T)
  fit.ridge = cv.glmnet(as.matrix(X.train.trim[, beta.snps]), y.train, alpha = 0, family = "gaussian")
  fit.df.ridge = data.frame(fit.ridge$cvm, fit.ridge$lambda, 1:length(fit.ridge$cvm))
  best.lam.ridge = fit.df.ridge[fit.df.ridge[,2] == fit.ridge$lambda.min, ]
  ret.ridge = as.data.frame(fit.ridge$glmnet.fit$beta[, best.lam.ridge[1,3]])
  ridge.cor = prettyNum(cor(as.matrix(X.train.trim[, beta.snps]) %*% ret.ridge[,1], y.train)^2, format = "fg", digits = 5)
  cat(paste0("Ridge model correlation squared is ", ridge.cor, "\n"))
  cat(paste0("Ridge model correlation squared is ", ridge.cor, "\n"), file = paste0(opt$o, ".log"), append = T)
  
  beta.values.ridge = ret.ridge[, 1]
  beta.snps.ridge = colnames(X.train.trim[, beta.snps])
  beta.snps.split = matrix(unlist(strsplit(beta.snps.ridge, split = ":")), nrow = length(beta.snps.ridge), ncol = 4, byrow = T)
  betas.ridge = cbind(beta.snps.split, beta.values.ridge)
  colnames(betas.ridge) = c("Chr", "Position", "Dosage_allele", "Other_allele", "Effect")
  write.table(betas.ridge, file = paste0(opt$o, ".betas.ridge.txt"), col.names = T, row.names = F, quote = F, sep = "\t")
} else {
  f = lm.summary$fstatistic
  p.val = prettyNum(1 - pf(f[1], f[2], f[3]), format = "fg", digits = 5)
  cat(paste0("P-value for linear model is ", p.val, "\n"))
  cat(paste0("Adjusted r squared for linear model is ", prettyNum(lm.summary$r.squared, format = "fg", digits = 5), "\n"))
  cat("\n\n************ Secondary Model Information ************\n", file = paste0(opt$o, ".log"), append = T)
  cat(paste0("P-value for linear model is ", p.val, "\n"), file = paste0(opt$o, ".log"), append = T)
  cat(paste0("Adjusted r squared for linear model is ", prettyNum(lm.summary$r.squared, format = "fg", digits = 5), "\n"), file = paste0(opt$o, ".log"), append = T)
  beta.snps.lm = gsub("`", "", rownames(lm.summary$coeff)[2:length(rownames(lm.summary$coeff))])
  beta.snps.split = matrix(unlist(strsplit(beta.snps.lm, split = ":")), nrow = length(beta.snps.lm), ncol = 4, byrow = T)
  betas.lm = cbind(rbind(c("Intercept", rep(NA, 3)), beta.snps.split), lm.summary$coeff[,1])
  colnames(betas.lm) = c("Chr", "Position", "Dosage_allele", "Other_allele", "Effect")
  write.table(betas.lm, file = paste0(opt$o, ".betas.lm.txt"), col.names = T, row.names = F, quote = F, sep = "\t")
}


## Stop timer and report total run time

script.time = proc.time() - time
cat(paste0("Total run time was ", script.time[3], " seconds\n\n"))
cat("\n\n************ Run Information ************\n", file = paste0(opt$o, ".log"), append = T)
cat(paste0("Total run time was ", script.time[3], " seconds\n"), file = paste0(opt$o, ".log"), append = T)