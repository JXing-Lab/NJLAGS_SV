# Figures and Tables

This directory contains three folders with the scripts used to generate the figures and tables for this project. Two of them contain the scripts for the CNV and gSV pipeline figures/tables respectively, while the third contains the scripts for generating the protein-protein interaction network used in Figure 5.

### CNV figures and tables

- Starting_Data_Scripts.py - Used to build Table 1 and Table 2
- CNV_Summarizer.py - Used to build Table 5
- CNV_Plots.py - Used to plot Figure S1(A)
- CNV_Enrichment.py - Used to build table for use in Figure S1(B)
- BoxPlotEnrichment.R - Used to plot Figure S1(B)

### gSV figures and tables
- 7_tables_and_figures.sh is used to generate the tables and figures
- GenerateGeneLevelSummaryTable_CNV.py - Used to build CNV input table for Table 4
- GenerateGeneLevelSummaryTable_MEISV.py - Used to build gSV input table for Table 4
- GenerateGeneLevelSummaryTable.py - Used to build Table 4
- GenerateSizeDistributionChartsByMEI.py - Used to plot Figure S2

### PPI Network

This folder contains the files and scripts used to generate Figure 5. Pathway data was drawn from ConsensusPathDB, and input genes from the results of the combined CNV/gSV pipeline.

- Network_gene.py - Used to build gene table referenced by Network_nodes_edges.py
- Network_nodes_edges.py - Used to build table defining network structure
- Network_plot.R - Used to plot PPI Network
