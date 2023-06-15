SRC_DIR=/lab01/Projects/Sammy_Projects/NJLAGS_SVs/src

# 1. Adding segregation data
python ${SRC_DIR}/GetSegregation.py

# 2. Adding Brain expression and IMPC data
python ${SRC_DIR}/GetBrainExpressionAndIMPC.py

# 3. Adding Sample AF data
python ${SRC_DIR}/GetSampAF.py

# 4. Adding eQTL data and Pop AF
python ${SRC_DIR}/GetSVeQTL.py
python ${SRC_DIR}/MergedeQTL.py
python ${SRC_DIR}/GeteQTL.py

# Removing mendel errors step is not needed, since I'm handling that in step 1.

# 5. Add NDD lists annotations
python ${SRC_DIR}/GetNDDs.py

# 6.Add biotype info
python ${SRC_DIR}/GetBiotype.py

