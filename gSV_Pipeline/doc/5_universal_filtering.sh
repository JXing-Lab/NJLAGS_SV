# Before splitting into seperate pipelines, do universal filtering on the fully annotated callset
SRC_DIR=/lab01/Projects/Sammy_Projects/NJLAGS_SVs/src

# 1. Inheritance Pattern Filtering
python ${SRC_DIR}/FilterInh.py

# 2. Benign Filtering
python ${SRC_DIR}/FilterBenign.py