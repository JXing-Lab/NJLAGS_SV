SRC_DIR=/lab01/Projects/Sammy_Projects/NJLAGS_SVs/src


# Part 1 - Exonic Pipeline
python ${SRC_DIR}/Pipelines/ExonicPipeline.py

# Part 2 - Intronic Pipelines
python ${SRC_DIR}/Pipelines/IntronicEQTLPipeline.py
python ${SRC_DIR}/Pipelines/IntronicAFPipeline.py

# Part 3 - Intergenic Pipeline
python ${SRC_DIR}/Pipelines/IntergenicPipeline.py

