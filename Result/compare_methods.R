setwd("/proj/yunligrp/users/chanhwa/pwas/whi_pwas_freeze9/method")
library(data.table)
library(rlang)
library(dplyr)
library(VennDiagram)
library(ggplot2)
library(ggrepel)
library(ggExtra)
library(cowplot)
library(tidyverse)
if (!require(devtools)) install.packages("devtools")
# devtools::install_github("gaospecial/ggVennDiagram")
library("ggVennDiagram")
if (!require(devtools)) install.packages("devtools")
# devtools::install_github("yanlinlin82/ggvenn")
library(ggvenn)
library(gtools)
# install.packages('Rcpp')
library(Rcpp)
library(gridExtra)

# Load rsq data
rsq = fread("correlations.txt")

# p-value thresholds
colnames(rsq)

# Number of failed model for each thresholds
colSums(apply(rsq, is.na, MARGIN = 2))


#' Draw Venn Diagram of successfully trained protein prediction models of two methods
#'
#' @param method1 A string. Column name of `rsq`.
#' @param method2 A string. Column name of `rsq`.
#' @param data A data frame. Data frame of rsq values of protein prediction models for various methods 
#' (cis-only, proposed method with various p-value thresholds).
#' @return Venn Diagram of the number of trained protein prediction models of two methods as well as
#' violin plot of rsq values of protein prediction models successfully trained by only one-method.
#' @examples 
#' 

venndiagram <- function(method1, method2, data){
  set1 <- data$protein[!is.na(data[[method1]])]
  set2 <- data$protein[!is.na(data[[method2]])]
  
  # Venn Diagram
  p1 <-  ggVennDiagram(x = list(set1, set2), category.names = c(method1 , method2), edge_size = 0.5) + 
    ggtitle("Number of successfully trained protein prediction model") +
    scale_fill_gradient(low="white",high = "red") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "none")
  
  # Summary statistics by category
  set1.only <- setdiff(set1, set2)
  set2.only <- setdiff(set2, set1)
  common <- intersect(set1, set2)
  
  summary.df <- as.data.frame(cbind(
    summary(data[[method1]][data$protein %in% set1.only]),
    summary(data[[method1]][data$protein %in% common]),
    summary(data[[method2]][data$protein %in% common]),
    summary(data[[method2]][data$protein %in% set2.only])
  ))
  
  colnames(summary.df) <- c(paste0(method1, " - only"), paste0(method1, " - common"), 
                            paste0(method2, " - common"), paste0(method2, " - only"))
  summary.df <- signif(summary.df, digits = 4)
  
  # rsq distribution by category
  idx.1 <- !is.na(data[[method1]])
  idx.2 <- !is.na(data[[method2]])
  
  data.cat = data %>% select(protein, {{method1}}, {{method2}}) %>% 
    mutate(method1.only = ifelse(idx.1 & !idx.2, .[[method1]], NA), 
           method1.common = ifelse(idx.1 & idx.2, .[[method1]], NA), 
           method2.common = ifelse(idx.1 & idx.2, .[[method2]], NA), 
           method2.only = ifelse(!idx.1 & idx.2, .[[method2]], NA)) %>%
    select(protein, method1.only, method1.common, method2.common, method2.only)
  
  colnames(data.cat) <- c("protein", paste0(method1, " - only"), paste0(method1, " - common"), 
                          paste0(method2, " - common"), paste0(method2, " - only"))
  
  data.melt = melt(data.cat, na.rm = T)
  
  p2 <- ggplot(data.melt, aes(x = variable, y = value)) + 
    geom_violin() +
    geom_boxplot(width = 0.05) +
    labs(x = "Method", y = "Testing Rsq") +
    ggtitle("Testing Rsq Violin Plots") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "none")
  
  plot_grid(p1, tableGrob(summary.df), p2, 
            nrow = 3,
            labels = c("", "Testing Rsq Summary statistics",""))
  
}

# Venn Diagrams for method pairs
methods = colnames(rsq)[-1]
methods_pair = combinations(n = length(methods), r = 2, v = methods)

pdf(file = "venndiagram.pdf",width = 10, height = 10)
for(i in 1:nrow(methods_pair)){
  print(methods_pair[i,])
  print(venndiagram(method1 = methods_pair[i,1], method2 = methods_pair[i,2], data = rsq))
}
dev.off()



#' Draw scatter plot of protein prediction model testing rsq of two methods
#'
#' @param method1 A string. Column name of `rsq`.
#' @param method2 A string. Column name of `rsq`.
#' @param data A data frame. Data frame of rsq values of protein prediction models for various methods 
#' (cis-only, proposed method with various p-value thresholds).
#' @param diff_threshold A number. Threshold that decide which method is dominant.
#' @return scatter plot of protein prediction model testing rsq of two methods.
#' @examples 
#' scatterplot(method1 = "cisonly", method2 = "proposed_0.001", data = rsq, threshold = 0.15)


scatterplot <- function(method1, method2, data, diff_thresh){
  data.comp = data %>% select(protein, {{method1}}, {{method2}})
  data.comp[["category"]] = "methods.sim"
  data.comp[["category"]][data.comp[[method1]] - data.comp[[method2]] > diff_thresh] = "method1.better"
  data.comp[["category"]][data.comp[[method2]] - data.comp[[method1]] > diff_thresh] = "method2.better"
  
  data.comp[["category"]][!is.na(data.comp[[method1]]) & is.na(data.comp[[method2]])] = "method1.only"
  method1.only.idx <- which(data.comp[["category"]] == "method1.only")
  method1.only.top.idx <- method1.only.idx[
      head(order(data.comp[[method1]][method1.only.idx],decreasing = T), min(5, length(method1.only.idx)))]
  data.comp[["category"]][method1.only.top.idx] = "method1.only.top"
  
  data.comp[["category"]][!is.na(data.comp[[method2]]) & is.na(data.comp[[method1]])] = "method2.only"
  method2.only.idx <- which(data.comp[["category"]] == "method2.only")
  method2.only.top.idx <- method2.only.idx[
    head(order(data.comp[[method2]][method2.only.idx],decreasing = T), min(5, length(method2.only.idx)))]
  data.comp[["category"]][method2.only.top.idx] = "method2.only.top"
  
  data.comp[["category"]][is.na(data.comp[[method2]]) & is.na(data.comp[[method1]])] = "no.methods"
  
  data.comp = data.comp %>% 
    mutate({{method1}} := ifelse( is.na(.[[method1]]) , -0.1, .[[method1]] ),
           {{method2}} := ifelse( is.na(.[[method2]]) , -0.1, .[[method2]] ))
  
  
  ggplot() +
    geom_vline(xintercept = 0, alpha = 0.7) +
    geom_hline(yintercept = 0) +
    xlim(-0.3,1) + ylim(-0.3,1)+
    geom_point(data = data.comp,
               mapping = aes(x = get(method1), y = get(method2)),
               color = ifelse(data.comp$category == "method1.better", "blue", 
                       ifelse(data.comp$category == "method2.better", "red", "black")),
               alpha = 0.5,
               size = 1) +
    geom_text_repel(data = data.comp, 
                    mapping = aes(x = get(method1), y=get(method2), 
                                  label = ifelse(category == "method1.better", protein, "")),
                    seed = 1,
                    force = 5,
                    nudge_x           = 1,
                    direction         = "y",
                    hjust             = 0,
                    segment.size      = 0.2,
                    segment.curvature = -0.1,
                    segment.angle     = 45,
                    segment.ncp       = 1,
                    segment.square    = TRUE,
                    segment.inflect   = FALSE,
                    xlim = c(NA, 1),
                    ylim = c(0, NA),
                    size = 3) +
    geom_text_repel(data = data.comp, 
                    mapping = aes(x = get(method1), y=get(method2), 
                                  label = ifelse(category == "method1.only.top", protein, "")),
                    seed = 1,
                    nudge_y           = -0.5,
                    direction         = "x",
                    hjust             = 0,
                    segment.size      = 0.2,
                    segment.curvature = -0.1,
                    segment.angle     = 45,
                    segment.ncp       = 1,
                    segment.square    = TRUE,
                    segment.inflect   = FALSE,
                    ylim = c(NA,-0.1),
                    size = 3,
                    angle = 90) +
    geom_text_repel(data = data.comp, 
                    mapping = aes(x = get(method1), y=get(method2), 
                                  label = ifelse(category == "method2.better", protein, "")),
                    seed = 1,
                    force = 5,
                    nudge_y           = 0.5,
                    direction         = "x",
                    vjust             = 0,
                    segment.size      = 0.2,
                    segment.curvature = -0.1,
                    segment.angle     = 45,
                    segment.ncp       = 1,
                    segment.square    = TRUE,
                    segment.inflect   = FALSE,
                    ylim = c(NA, 1),
                    xlim = c(0, NA),
                    size = 3,
                    angle = 90) +
    geom_text_repel(data = data.comp, 
                    mapping = aes(x = get(method1), y=get(method2), 
                                  label = ifelse(category == "method2.only.top", protein, "")),
                    seed = 1,
                    nudge_x           = -0.5,
                    direction         = "y",
                    hjust             = 0,
                    segment.size      = 0.2,
                    segment.curvature = -0.1,
                    segment.angle     = 45,
                    segment.ncp       = 1,
                    segment.square    = TRUE,
                    segment.inflect   = FALSE,
                    xlim = c(NA,-0.1),
                    size = 3) +
    geom_abline(slope = 1, intercept = 0, color="red") +
    geom_abline(slope = 1, intercept = diff_thresh, color="red", lty = "dashed") +
    geom_abline(slope = 1, intercept = -diff_thresh, color="blue", lty = "dashed") +
    labs(x = paste0(method1, "_R2"), y = paste0(method2, "_R2")) +
    ggtitle(paste0("Scatter Plot of Testing Rsqs of protein prediction models from two methods\n",
                   "(",method1," vs ",method2,": Difference threshold = ",diff_thresh, ")")) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "none")
}






# #----- Older Version -----#
# scatterplot <- function(method1, method2, data, diff_thresh){
#   
#   data.sub = data %>% select(protein, all_of(method1), all_of(method2))
#   data.comp = data.sub[complete.cases(data.sub),]
#   data.comp = data.comp %>% 
#     mutate(dominate = ifelse(get(method1) - get(method2) > diff_thresh, method1, 
#                              ifelse(get(method2) - get(method1) > diff_thresh, method2, "none")))
#   
#   
#   ggplot(data = data.comp, 
#          aes(x = get(method1), y = get(method2),
#              label = ifelse(dominate %in% c(method1, method2), protein, ""))) +
#     geom_text_repel(seed = 1,
#                     nudge_x           = 0.25,
#                     direction         = "y",
#                     hjust             = 0,
#                     segment.size      = 0.2,
#                     segment.curvature = -0.1,
#                     segment.angle     = 45,
#                     segment.ncp       = 1,
#                     segment.square    = TRUE,
#                     segment.inflect   = FALSE,
#                     xlim = c(0.5,NA),
#                     size = 3) +
#     geom_point(color = ifelse(data.comp$dominate == method1, "blue", 
#                               ifelse(data.comp$dominate == method2, "red", "black"))) +
#     geom_abline(slope = 1, intercept = 0, color="red") +
#     geom_abline(slope = 1, intercept = diff_thresh, color="red", lty = "dashed") +
#     geom_abline(slope = 1, intercept = -diff_thresh, color="blue", lty = "dashed") +
#     xlim(-0.05,0.75) + ylim(-0.05,0.75)+
#     labs(x = paste0(method1, "_R2"), y = paste0(method2, "_R2")) +
#     ggtitle(paste0("Scatter Plot of Testing Rsqs of protein prediction models from two methods\n",
#                    "(",method1," vs ",method2,": Difference threshold = ",diff_thresh, ")")) +
#     theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "none")
# }
#----- Older Version -----#



#--- Scatter plots for method pairs ---#
methods = colnames(rsq)[-1]
methods_pair = combinations(n = length(methods), r = 2, v = methods)

pdf(file = "scatterplots.pdf")
for(i in 1:nrow(methods_pair)){
  print(methods_pair[i,])
  print(scatterplot(method1 = methods_pair[i,1], method2 = methods_pair[i,2], data = rsq, diff_thresh = 0.15))
}
dev.off()


#' Violin plots and Summary statistics of protein prediction model testing rsq available for all methods
#'
#' @param data A data frame. Data frame of rsq values of protein prediction models for various methods 
#' (cis-only, proposed method with various p-value thresholds).
#' @param threshold A number. Threshold that rsq smaller than this value will be filtered out.
#' @return Violin plots and Summary statistics of protein prediction model testing rsq of methods.

violinplot <- function(data, threshold){
  
  data.filter[data.filter < threshold] = NA
  data.filter = data.filter[complete.cases(data.filter)]
  
  # Summary statistics by methods
  # summary.df = c()
  methods = colnames(data.filter)[-1]
  # for(method in methods) summary.df = cbind(summary.df, summary(data.filter[[method]]))
  summary.df = sapply(methods, function(method) summary(data.filter[[method]]))
  summary.df <- as.data.frame(signif(summary.df, digits = 4))
  colnames(summary.df) <- methods
  
  # Violin plots by methods
  data.filter.melt = melt(data.filter)
  p1 <- ggplot(data.filter.melt, aes(x = variable, y = value)) + 
    geom_violin() +
    geom_boxplot(width = 0.05) +
    labs(x = "Method", y = "Testing Rsq") +
    ggtitle(paste0("Violin Plots of Testing Rsqs of protein prediction models from methods\n", "(Rsq threshold = ",threshold, ")")) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "none")
  
  # Violon plots and Summary statistics
  plot_grid(p1, tableGrob(summary.df), 
            nrow = 2)
  
}

#--- Rsq Violin plots ---#
pdf(file = "violinplots.pdf",width = 13, height = 10)
for(threshold in c(0.001, 0.01, 0.05, 0.1, 0.2, 0.5)){
  print(threshold)
  print(violinplot(data = rsq, threshold = threshold))
}
dev.off()



















#' Summary statistics for proteins that are in the middle part of the venn diagram for every methods
#'
#' @param data A data frame. Data frame of rsq values of protein prediction models for various methods 
#' (cis-only, proposed method with various p-value thresholds).
#' @param threshold A number. Threshold that rsq smaller than this value will be filtered out.
#' @return Violin plots of protein prediction model testing rsq of methods.

violinplot <- function(data, threshold){
  
  data.filter = data
  
  data.filter %>% select()
  
  for(method in methods){
    data.filter[[method]][data.filter[[method]] < threshold] = NA
  }
  
  data.filter.melt = melt(data.filter)
  
  ggplot(data.filter.melt, aes(x = variable, y = value)) + 
    geom_violin() +
    geom_boxplot(width = 0.05) +
    labs(x = "Method", y = "Testing Rsq") +
    ggtitle(paste0("Violin Plots of Testing Rsqs of protein prediction models (available by all methods)\n", 
                   "(Rsq threshold = ",threshold, ")")) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "none")
  
}

#--- Rsq Violin plots ---#
pdf(file = "violinplots.pdf",width = 10, height = 6)
for(threshold in c(0.001, 0.01, 0.05, 0.1, 0.2, 0.5)){
  print(threshold)
  print(violinplot(data = rsq, threshold = threshold))
}
dev.off()

