# Libraries
library(ggplot2)
library(dplyr)
library(tidyverse)
library(hrbrthemes)
library(viridis)
library(Cairo)

# Read in the data
setwd("/home/rohan/CNV_Project/2023/Submission_Pipeline/Temp/Results/Plots/Enrichment")
dfu <- read.csv("unaffected_CDS_lengths.tsv", sep = "\t")
dfp_ASD <- read.csv("ASD_proband_CDS_lengths.tsv", sep = "\t") 
dfp_LI <- read.csv("LI_proband_CDS_lengths.tsv", sep = "\t") 
dfp_RI <- read.csv("RI_proband_CDS_lengths.tsv", sep = "\t") 

# calculate medians
unaff_median <- median(dfu$CDS_impact)
ASD_proband_median <- median(dfp_ASD$CDS_impact)
LI_proband_median <- median(dfp_LI$CDS_impact)
RI_proband_median <- median(dfp_RI$CDS_impact)

# print out medians
print(unaff_median)
print(ASD_proband_median)
print(LI_proband_median)
print(RI_proband_median)

# calculate means
unaff_mean <- mean(dfu$CDS_impact)
ASD_proband_mean <- mean(dfp_ASD$CDS_impact)
LI_proband_mean <- mean(dfp_LI$CDS_impact)
RI_proband_mean <- mean(dfp_RI$CDS_impact)

# print out means
print(unaff_median)
print(ASD_proband_median)
print(LI_proband_median)
print(RI_proband_median)

# set lists
unaff_CDS_lengths <- list(dfu$CDS_impact)
ASD_proband_CDS_lengths <- list(dfp_ASD$CDS_impact)
LI_proband_CDS_lengths <- list(dfp_LI$CDS_impact)
RI_proband_CDS_lengths <- list(dfp_RI$CDS_impact)

# prep data for graph
unaff = data.frame(group = "unaffected", value = unaff_CDS_lengths)
colnames(unaff) <- c("group", "value")
ASD_proband = data.frame(group = "ASD patients", value = ASD_proband_CDS_lengths)
colnames(ASD_proband) <- c("group", "value")
LI_proband = data.frame(group = "LI patients", value = LI_proband_CDS_lengths)
colnames(LI_proband) <- c("group", "value")
RI_proband = data.frame(group = "RI patients", value = RI_proband_CDS_lengths)
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

ggsave("CNV_burden_CDS_box.png")
ggsave("CNV_burden_CDS_box.pdf", device = cairo_pdf)

# plot data in violin
ggplot(plot.data, aes(x = group, y = value, fill = group)) +
  geom_violin() +
  ylab("CNVs") +
  theme(text=element_text(size=20,  family="Times New Roman")) +
  scale_fill_manual(values=c("#ff6464", "#1b92e8", "#e8de1b", "#1be826"))
#+ geom_jitter(color="black", size=0.5, alpha=0.9)

ggsave("CNV_burden_CDS_violin.png")
ggsave("CNV_burden_CDS_violin.pdf", device = cairo_pdf)

# run chi-squared tests
#m <- cbind(table(dfp_ASD$CDS_impact), table(dfu$CDS_impact))
#chisq.test(m)
#m <- cbind(table(dfp_LI$CDS_impact), table(dfu$CDS_impact))
#chisq.test(m)
#m <- cbind(table(dfp_RI$CDS_impact), table(dfu$CDS_impact))
#chisq.test(m)

# run Mann-Whitney-Wilcoxon tests
unaff_CDS_lengths <- as.numeric(unlist(unaff_CDS_lengths))
ASD_proband_CDS_lengths <- as.numeric(unlist(ASD_proband_CDS_lengths))
ASD_res = wilcox.test(unaff_CDS_lengths, ASD_proband_CDS_lengths, alternative = "two.sided")
ASD_res
LI_proband_CDS_lengths <- as.numeric(unlist(LI_proband_CDS_lengths))
LI_res = wilcox.test(unaff_CDS_lengths, LI_proband_CDS_lengths, alternative = "two.sided")
LI_res
RI_proband_CDS_lengths <- as.numeric(unlist(RI_proband_CDS_lengths))
RI_res = wilcox.test(unaff_CDS_lengths, RI_proband_CDS_lengths, alternative = "two.sided")
RI_res
