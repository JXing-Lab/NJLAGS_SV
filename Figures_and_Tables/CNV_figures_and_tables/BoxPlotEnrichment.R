# Libraries
library(ggplot2)
library(dplyr)
library(tidyverse)
library(hrbrthemes)
library(viridis)
library(Cairo)

# Read in the data
setwd("/home/rohan/CNV_Project/2023/Submission_Pipeline/Temp/Results/Plots/Enrichment")
dfu <- read.csv("unaffected_counts.tsv", sep = "\t")
dfp_ASD <- read.csv("ASD_proband_counts.tsv", sep = "\t") 
dfp_LI <- read.csv("LI_proband_counts.tsv", sep = "\t") 
dfp_RI <- read.csv("RI_proband_counts.tsv", sep = "\t") 

# calculate medians
unaff_median <- median(dfu$CN_sum)
ASD_proband_median <- median(dfp_ASD$CN_sum)
LI_proband_median <- median(dfp_LI$CN_sum)
RI_proband_median <- median(dfp_RI$CN_sum)

# print out medians
print(unaff_median)
print(ASD_proband_median)
print(LI_proband_median)
print(RI_proband_median)

# calculate means
unaff_mean <- mean(dfu$CN_sum)
ASD_proband_mean <- mean(dfp_ASD$CN_sum)
LI_proband_mean <- mean(dfp_LI$CN_sum)
RI_proband_mean <- mean(dfp_RI$CN_sum)

# print out means
print(unaff_median)
print(ASD_proband_median)
print(LI_proband_median)
print(RI_proband_median)

# set lists
unaff_counts <- list(dfu$CN_sum)
ASD_proband_counts <- list(dfp_ASD$CN_sum)
LI_proband_counts <- list(dfp_LI$CN_sum)
RI_proband_counts <- list(dfp_RI$CN_sum)

# prep data for graph
unaff = data.frame(group = "unaffected", value = unaff_counts)
colnames(unaff) <- c("group", "value")
ASD_proband = data.frame(group = "ASD patients", value = ASD_proband_counts)
colnames(ASD_proband) <- c("group", "value")
LI_proband = data.frame(group = "LI patients", value = LI_proband_counts)
colnames(LI_proband) <- c("group", "value")
RI_proband = data.frame(group = "RI patients", value = RI_proband_counts)
colnames(RI_proband) <- c("group", "value")

plot.data = rbind(unaff, ASD_proband, LI_proband, RI_proband)

# plot data in boxplot
ggplot(plot.data, aes(x = group, y = value, fill = group)) +
  geom_boxplot(varwidth = TRUE) +
  ylab("CNVs") +
  theme(text=element_text(size=20,  family="Times New Roman"), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  scale_fill_manual(values=c("#ff6464", "#1b92e8", "#e8de1b", "#1be826")) +  
  geom_jitter(color="black", size=0.5, alpha=0.9)

ggsave("CNV_burden_box.png")
ggsave("CNV_burden_box.pdf", device = cairo_pdf)

# plot data in violin
ggplot(plot.data, aes(x = group, y = value, fill = group)) +
  geom_violin() +
  ylab("CNVs") +
  theme(text=element_text(size=20,  family="Times New Roman")) +
  scale_fill_manual(values=c("#ff6464", "#1b92e8", "#e8de1b", "#1be826"))
  #+ geom_jitter(color="black", size=0.5, alpha=0.9)

ggsave("CNV_burden_violin.png")
ggsave("CNV_burden_violin.pdf", device = cairo_pdf)

# run chi-squared tests
m <- cbind(table(dfp_ASD$CN_sum), table(dfu$CN_sum))
chisq.test(m)
m <- cbind(table(dfp_LI$CN_sum), table(dfu$CN_sum))
chisq.test(m)
m <- cbind(table(dfp_RI$CN_sum), table(dfu$CN_sum))
chisq.test(m)

# run Mann-Whitney-Wilcoxon tests
unaff_counts <- as.numeric(unlist(unaff_counts))
ASD_proband_counts <- as.numeric(unlist(ASD_proband_counts))
ASD_res = wilcox.test(unaff_counts, ASD_proband_counts, alternative = "two.sided")
ASD_res
LI_proband_counts <- as.numeric(unlist(LI_proband_counts))
LI_res = wilcox.test(unaff_counts, LI_proband_counts, alternative = "two.sided")
LI_res
RI_proband_counts <- as.numeric(unlist(RI_proband_counts))
RI_res = wilcox.test(unaff_counts, RI_proband_counts, alternative = "two.sided")
RI_res
