#!/bin/bash
#$ -V
#$ -S /bin/bash

source ./pipeline.cfg

# get input
input_dir="${work_dir}/bam/"
input=""
while read file_name; do
  input="${input}-I ${input_dir}${file_name}.recal.bam "
done < $1

for i in {1..12}; do
  input="${input}-I ${input_dir}${i}.recal.bam "
done

# log function
time_start=$(date +%s)
time_now=${time_start}
time_end=${time_start}
LOG_FILE="${work_dir}/logs/${file_name}-${time_start}.log"
LOG() {
  time_now=$(date +%s)
  echo "[${file_name}][$(date "+%F %T")][time: $((${time_now} - ${time_end})) sec] $*" >> ${LOG_FILE}
  time_end=${time_now}
}


LOG "Start processing sample ${file_name}..."

# 8
echo "java -Xmx15g -XX:ParallelGCThreads=8 -Djava.io.tmpdir=${tmp} -jar ${gatk} -T UnifiedGenotyper -nt 8 -R ${ref} -I ${input} -o "${work_dir}/112.vcf" -glm BOTH -dcov 200 -stand_call_conf 30 -stand_emit_conf 30 --dbsnp ${dbsnp}"
java -Xmx15g -XX:ParallelGCThreads=8 -Djava.io.tmpdir=${tmp} -jar ${gatk} -T UnifiedGenotyper -nt 8 -R ${ref} -I ${input} -o "${work_dir}/112.vcf" -glm BOTH -dcov 200 -stand_call_conf 30 -stand_emit_conf 30 --dbsnp ${dbsnp}

LOG "Done UnifiedGenotyper, total elasped time: $((${time_now} - ${time_start})) sec."
