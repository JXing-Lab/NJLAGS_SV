
# CNV Pipeline

### Input files

Data that has to be input into the pipeline from external sources, rather than being generated by the pipeline itself.

Builder phase:

- starting_data_files: microarray files from Mahaganapathy
- NJLAGS_CNV.ped: pedigree file of NJLAGS cohort
- ASD_only.ped: pedigree file for ASD_only phenotype
- ASD_LI.ped: pedigree file for ASD_only phenotype
- ASD_RI.ped: pedigree file for ASD_only phenotype
- human_g1k_v37.fasta: used for conversion to VCF

Prioritization phase:

- NDD_genes.txt: neurodevelopmental disorder genes, from SFARI
- dispensable_genes_and_muc: genes to exclude, from Rausell et al
- GTEx_2019_12_12.txt: gene expression file from GTEx
- tpm1.txt/tpm2.txt/tpm3.txt: gene expression files from other sources
- syndromes.bed: known ASD related syndromes, from SFARI

### Instructions

The CNV_Master.py script will run the entire pipeline from start to finish. It does so across three phases, calling the corresponding Python scripts as it does so:

- CNV_Builder.py
    - Merge CNV files across batches and callers
    - Convert to VCF format for use in Annotation
- CNV_Annotator.py
    - QC filtration to remove outliers
    - Running AnnotSV annotation on VCF
    - Running GEMINI-based Python script “Geminesque.py” on annotated VCF
    - Filtering on segregation analysis
- CNV_Prioritizer.py
    - Prioritizing candidate variants on dispensability, coding region overlap, gene expression, internal cohort frequency
    - Calling StrVCTVRE to predict pathogenicity
    - Calling SvAnna to predict pathogenicity

Additionally, the CNV_Summarizer.py script creates summary tables for analysis.

### Running CNV_Master.py

CNV_Master.py already has all the commands for the constituent scripts. All that needs to be specified is the project folder where the pipeline is being run and the proband phenotype of interest.

```python
# DECLARE GLOBAL VARIABLES
projects_folder = "PROJECT FOLDER FILEPATH"
phenotype = "PHENOTYPE" # set phenotype
print("Projects folder located at: " + projects_folder)
```