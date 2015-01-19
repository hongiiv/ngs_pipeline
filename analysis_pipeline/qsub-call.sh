#!/bin/bash
# $1 - file_name-sample_id mapping file path

source ./pipeline.cfg

while read file_name
qsub -pe make 8 -cwd -j y ./pipeline_unified_genotyper.sh $1
