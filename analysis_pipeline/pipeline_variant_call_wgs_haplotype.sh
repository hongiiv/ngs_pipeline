#!/bin/bash
#$ -V
#$ -S /bin/bash

# Haplotype Caller script for WGS sample
# runs on 8vCore, 64GB RAM machine

source ./pipeline.cfg

# arguments
sample_id=$1

# log function
time_start=$(date +%s)
time_now=${time_start}
time_end=${time_start}
LOG_FILE="./logs/kr37wgs.log"
LOG() {
  time_now=$(date +%s)
  echo "[KR37wgs][$(date "+%F %T")][time: $((${time_now} - ${time_end})) sec] $*" >> ${LOG_FILE}
  time_end=${time_now}
}

# variables
bam_dir="${work_dir}results/kr37/bam/"
out_dir="${work_dir}"

mkdir ${tmp}

LOG "Start processing KR37wgs..."

## Haplotype Caller for each sample
echo "java -Xmx60g -XX:ParallelGCThreads=8 -Djava.io.tmpdir=${tmp} -jar ${gatk} -T HaplotypeCaller -R ${ref} -I ${bam_dir}${sample_id}.recal.bam --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 -stand_call_conf 30.0 -stand_emit_conf 10.0 -o ${out_dir}${sample_id}.hc.vcf"
java -Xmx60g -XX:ParallelGCThreads=8 -Djava.io.tmpdir=${tmp} -jar ${gatk} -T HaplotypeCaller -R ${ref} -I ${bam_dir}${sample_id}.recal.bam --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 -stand_call_conf 30.0 -stand_emit_conf 10.0 -o ${out_dir}${sample_id}.hc.vcf
