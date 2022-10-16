library(tidyverse)

data <- read.csv("olink_data_2020SEP10.csv", header = T)
proteins <- colnames(data)[-1:-7]
proteins <- proteins[-(1+(0:5)*93)]
olink.proteins <- data.frame(protein = unlist(strsplit(x = proteins, split = "_"))[1:552*2],
               assay = unlist(strsplit(x = proteins, split = "_"))[-1+1:552*2])

protein_info <- read.csv("olink_proteins.csv") %>% mutate(protein = gsub("[ -]",".",protein))     

result <- merge(x = olink.proteins, y = protein_info, by = "protein", all.x = T)
(result %>% filter(is.na(assay.y)))$protein

right_result <- merge(x = olink.proteins, y = protein_info, by = "protein", all.y = T)
right_result %>% filter(is.na(assay.x))

olink_protein_name_nomatch <- c("ANGPT1", "BOC", "ICAM.5", "IL12", "IL13", "IL18", "IL18", "IL2", "IL33", "IL4", "IL7", "IL8", "MMP12", "MMP7", "PRSS8", "TXNDC15", "VEGFA", "X4E.BP1")

protein_info_name_nomatch <- c("ANG.1","Protein.BOC","ICAM5","IL.12","IL.13","IL.18","IL.18","IL.2","IL.33","IL.4","IL.7","IL.8","MMP.12","MMP.7","PRSS8.","TXD15","VEGF.A","4E.BP1")

for(i in 1:length(protein_info$protein)){
  match_idx <- match(protein_info$protein[i], protein_info_name_nomatch)
  if(!is.na(match_idx)) protein_info$protein[i] <- olink_protein_name_nomatch[match_idx]
}

result <- merge(x = olink.proteins, y = protein_info, by = c("protein", "assay"))

begin <- as.vector(t(result %>% filter(assay == "cardio2") %>% select(begin)))
end <- as.vector(t(result %>% filter(assay == "cardio2") %>% select(end)))

result[which(result$assay == "cardio2"), ]$begin <- end
result[which(result$assay == "cardio2"), ]$end <- begin

output <- result %>% mutate(protein = paste(assay, protein,sep = "_")) %>% select(protein, chr, begin, end)

write.table(output,"proteins_info.txt", sep = "\t", row.names = F, quote = F, col.names = F)


###############################################################################################

cis_corr <- read.delim("correlations.txt", sep = " ", header = F)
colnames(cis_corr) <- c("protein", "corr")
no_corr <- cis_corr %>% filter(is.na(corr)) # %>% select(protein)
no_corr_output <- merge(x = output, y = no_corr, by = "protein") %>% select(protein, chr, begin, end)
write.table(output,"no_corr_proteins_info.txt", sep = "\t", row.names = F, quote = F, col.names = F)

no_model_conv <- as.vector(t(read.delim("no_corr_reason.txt", header = F)))
redo <- no_corr_output %>% filter(!(protein %in% no_model_conv))
write.table(redo,"redo_proteins.txt", sep = "\t", row.names = F, quote = F, col.names = F)
