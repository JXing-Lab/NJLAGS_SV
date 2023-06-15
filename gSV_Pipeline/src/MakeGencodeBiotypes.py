from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation
from Bio.SeqFeature import SeqFeature
import pandas as pd

annotation_file = "/lab01/DataSets/hg19_annotation/gencode_v14/gencode.v14.annotation.gtf"  # Replace with the path to your annotation file

data = []
with open(annotation_file, "r") as file:
    for line in file:
        if line.startswith("#"):
            continue
        fields = line.strip().split("\t")
        attributes = dict(item.strip().split(" ") for item in fields[8].split(";") if item.strip())
        gene_type = attributes.get("gene_type", "")
        gene_name = attributes.get("gene_name", "")
        if gene_type == '"protein_coding"':
            data.append([gene_type, gene_name])

df = pd.DataFrame(data, columns=["Gene Type", "Gene Name"])
df=df.drop_duplicates()
df["Gene Type"]=df["Gene Type"].str.strip('""')
df["Gene Name"]=df["Gene Name"].str.strip('""')

df.to_csv("/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/gencode_biotypes.csv", index=None)
