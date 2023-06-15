SRC_DIR=/lab01/Projects/Sammy_Projects/NJLAGS_SVs/src

# Counts for each unique gene and fam for each phenotype
python ${SRC_DIR}/tables_and_figures/GenerateCountsByPhenotypeTable.py

# Phenotype by sex
python ${SRC_DIR}/tables_and_figures/GeneratePhenotypeCountsBySex.py

# Size distrib chart and figure
python ${SRC_DIR}/tables_and_figures/GenerateSizeDistributionChartsByMEI.py

# MEI/SV Variant-level summary table
python ${SRC_DIR}/tables_and_figures/GenerateVariantLevelSummaryTable.py

# CNV/MEI/SV Gene-level summary table
python ${SRC_DIR}/tables_and_figures/GenerateGeneLevelSummaryTable_MEISV.py
python ${SRC_DIR}/tables_and_figures/GenerateGeneLevelSummaryTable_CNV.py
python ${SRC_DIR}/tables_and_figures/GenerateGeneLevelSummaryTable.py

# Missing rates for each individual from MEI files
python ${SRC_DIR}/tables_and_figures/GenerateMEIMissingRateByIndividualTable.py

# Inheritance pattern counts by phenotype
python ${SRC_DIR}/tables_and_figures/GenerateInheritanceCountsByPhenotypeTable.py

# Counts of various dataframes to assist in making the big flowchart
python ${SRC_DIR}/tables_and_figures/GenerateDataframeCountsTable.py

# Make SV Burden box and whisker plots plus stats
python ${SRC_DIR}/tables_and_figures/GenerateSVBurdenFigure.py