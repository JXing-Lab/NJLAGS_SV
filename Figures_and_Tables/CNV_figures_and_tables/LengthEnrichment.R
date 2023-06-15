# Libraries
library(ggplot2)
library(dplyr)
library(tidyverse)
library(hrbrthemes)
library(viridis)
library(Cairo)

# Read in the data
setwd("/home/rohan/CNV_Project/2023/Submission_Pipeline/Temp/Results/Plots/Enrichment")
dfu <- read.csv("unaffected_lengths.tsv", sep = "\t")
dfp_ASD <- read.csv("ASD_proband_lengths.tsv", sep = "\t") 
dfp_LI <- read.csv("LI_proband_lengths.tsv", sep = "\t") 
dfp_RI <- read.csv("RI_proband_lengths.tsv", sep = "\t") 

# calculate medians
unaff_median <- median(dfu$length_impact_1)
ASD_proband_median <- median(dfp_ASD$length_impact_1)
LI_proband_median <- median(dfp_LI$length_impact_1)
RI_proband_median <- median(dfp_RI$length_impact_1)

# print out medians
print(unaff_median)
print(ASD_proband_median)
print(LI_proband_median)
print(RI_proband_median)

# calculate means
unaff_mean <- mean(dfu$length_impact_1)
ASD_proband_mean <- mean(dfp_ASD$length_impact_1)
LI_proband_mean <- mean(dfp_LI$length_impact_1)
RI_proband_mean <- mean(dfp_RI$length_impact_1)

# print out means
print(unaff_median)
print(ASD_proband_median)
print(LI_proband_median)
print(RI_proband_median)

# set lists
unaff_lengths <- list(dfu$length_impact_1)
ASD_proband_lengths <- list(dfp_ASD$length_impact_1)
LI_proband_lengths <- list(dfp_LI$length_impact_1)
RI_proband_lengths <- list(dfp_RI$length_impact_1)

# prep data for graph
unaff = data.frame(group = "unaffected", value = unaff_lengths)
colnames(unaff) <- c("group", "value")
ASD_proband = data.frame(group = "ASD patients", value = ASD_proband_lengths)
colnames(ASD_proband) <- c("group", "value")
LI_proband = data.frame(group = "LI patients", value = LI_proband_lengths)
colnames(LI_proband) <- c("group", "value")
RI_proband = data.frame(group = "RI patients", value = RI_proband_lengths)
colnames(RI_proband) <- c("group", "value")

plot.data = rbind(unaff, ASD_proband, LI_proband, RI_proband)

# plot data in boxplot
ggplot(plot.data, outlier.shape = NA, aes(x = group, y = value, fill = group)) +
  geom_boxplot(varwidth = TRUE) +
  scale_y_log10() +
  ylab("CNVs") +
  theme(text=element_text(size=20,  family="Times New Roman"), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  scale_fill_manual(values=c("#ff6464", "#1b92e8", "#e8de1b", "#1be826")) +  
  geom_jitter(color="black", size=0.5, alpha=0.9)

ggsave("CNV_burden_length_box.png")
ggsave("CNV_burden_length_box.pdf", device = cairo_pdf)

# plot data in violin
ggplot(plot.data, aes(x = group, y = value, fill = group)) +
  geom_violin() +
  ylab("CNVs") +
  theme(text=element_text(size=20,  family="Times New Roman")) +
  scale_fill_manual(values=c("#ff6464", "#1b92e8", "#e8de1b", "#1be826"))
#+ geom_jitter(color="black", size=0.5, alpha=0.9)

ggsave("CNV_burden_length_violin.png")
ggsave("CNV_burden_length_violin.pdf", device = cairo_pdf)

# run chi-squared tests
m <- cbind(table(dfp_ASD$length_impact_1), table(dfu$length_impact_1))
chisq.test(m)
m <- cbind(table(dfp_LI$length_impact_1), table(dfu$length_impact_1))
chisq.test(m)
m <- cbind(table(dfp_RI$length_impact_1), table(dfu$length_impact_1))
chisq.test(m)

# run Mann-Whitney-Wilcoxon tests
unaff_lengths <- as.numeric(unlist(unaff_lengths))
ASD_proband_lengths <- as.numeric(unlist(ASD_proband_lengths))
ASD_res = wilcox.test(unaff_lengths, ASD_proband_lengths, alternative = "two.sided")
ASD_res
LI_proband_lengths <- as.numeric(unlist(LI_proband_lengths))
LI_res = wilcox.test(unaff_lengths, LI_proband_lengths, alternative = "two.sided")
LI_res
RI_proband_lengths <- as.numeric(unlist(RI_proband_lengths))
RI_res = wilcox.test(unaff_lengths, RI_proband_lengths, alternative = "two.sided")
RI_res
