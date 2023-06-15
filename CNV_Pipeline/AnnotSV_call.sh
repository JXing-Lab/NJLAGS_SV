#!/bin/sh

export ANNOTSV=/lab01/Tools/AnnotSV
$ANNOTSV/bin/AnnotSV -SVinputFile "../Data/CNV3.vcf" -outputFile "../Data/CNV3_anno"
