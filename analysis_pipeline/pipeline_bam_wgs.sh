#!/bin/bash
#$ -V
#$ -S /bin/bash

# BAM finishing script for WGS sample
# runs on 8vCore, 16GB RAM machine

source ./pipeline.cfg

# arguments
category=$1
file_name=$2
sample_id=$3

# log function
time_start=$(date +%s)
time_now=${time_start}
time_end=${time_start}
LOG_FILE="./logs/${sample_id}.log"
LOG() {
  time_now=$(date +%s)
  echo "[${sample_id}][$(date "+%F %T")][time: $((${time_now} - ${time_end})) sec] $*" >> ${LOG_FILE}
  time_end=${time_now}
}

# functions 
check_error() {
if [ $1 -ne 0 ]; then
  echo "ERROR: pipeline.sh"
  echo "ERROR CODE: $1"
fi
}

check_line_size_bam() {
  fastqR1linesize=`zcat $1 | wc -l | cut -d' ' -f1`
  check_error $?

  fastqhalfsize=`expr ${fastqR1linesize} / 2`
  check_error $?

  bamlinesize=`${samtools} view $2 | wc -l ${FASTQ} | cut -d' ' -f1`
  check_error $?

  LOG "$1 half line size: ${fastqhalfsize}"
  LOG "$2 line size: ${bamlinesize}"
  if [ ${fastqhalfsize} -ne ${bamlinesize} ]; then
    LOG "check_line_size_bam failed!"
  fi
}

# variables
fastq_r1="${work_dir}fastq/disk6/${file_name}_R1.fastq.gz"
fastq_r2="${work_dir}fastq/disk6/${file_name}_R2.fastq.gz"
output="${work_dir}results/${category}/bam/${sample_id}"
stat="${work_dir}results/${category}/stat/${sample_id}"
#fastqc="${work_dir}results/${category}/fastqc/"
sam_path="/SNUH/SNUH_2/results/${category}/bam/${sample_id}"

LOG "Start processing sample ${sample_id}..."

# 0 Stat
(/SNUH/app/FastQC/fastqc -t 4 -o ${fastqc} -f fastq ${fastq_r1})&
(/SNUH/app/FastQC/fastqc -t 4 -o ${fastqc} -f fastq ${fastq_r2})&
wait
LOG "#0 FastQC done."

mv ${fastqc}${file_name}"_R1_fastqc.zip" ${fastqc}${sample_id}"_R1_fastqc.zip"
mv ${fastqc}${file_name}"_R2_fastqc.zip" ${fastqc}${sample_id}"_R2_fastqc.zip"

# 1
echo "#1"
echo "${bwa} mem -M -t 8 ${ref} ${fastq_r1} ${fastq_r2} > ${output}.sam"
${bwa} mem -M -t 8 ${ref} ${fastq_r1} ${fastq_r2} > ${output}.sam
LOG "#1 BWA done."

# 2
echo "#2"
LOG "#2 AddOrReplaceReadGroups start."
echo "java -Xmx15g -XX:ParallelGCThreads=8 -Djava.io.tmpdir=${tmp} -jar ${picard_folder}AddOrReplaceReadGroups.jar MAX_RECORDS_IN_RAM=350000 INPUT=${sam_path}.sam OUTPUT=${output}.sorted.bam SORT_ORDER=coordinate RGID=${sample_id} RGLB=${sample_id} RGPL=illumina RGPU=SureSelectAllExon RGSM=${sample_id} VALIDATION_STRINGENCY=LENIENT"
java -Xmx15g -XX:ParallelGCThreads=8 -Djava.io.tmpdir=${tmp} -jar ${picard_folder}AddOrReplaceReadGroups.jar MAX_RECORDS_IN_RAM=350000 INPUT=${sam_path}.sam OUTPUT=${output}.sorted.bam SORT_ORDER=coordinate RGID=${sample_id} RGLB=${sample_id} RGPL=illumina RGPU=SureSelectAllExon RGSM=${sample_id} VALIDATION_STRINGENCY=LENIENT
LOG "#2 AddOrReplaceReadGroups done."
check_line_size_bam ${fastq_r1} ${output}.sorted.bam

# 3
echo "#3"
LOG "#3 MarkDuplicates start."
echo "java -Xmx15g -XX:ParallelGCThreads=8 -Djava.io.tmpdir=${tmp} -jar ${picard_folder}MarkDuplicates.jar MAX_RECORDS_IN_RAM=350000 INPUT=${output}.sorted.bam OUTPUT=${output}.sorted.dp.bam METRICS_FILE=${output}.sorted.dp.metrix ASSUME_SORTED=true REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true"
java -Xmx15g -XX:ParallelGCThreads=8 -Djava.io.tmpdir=${tmp} -jar ${picard_folder}MarkDuplicates.jar MAX_RECORDS_IN_RAM=350000 INPUT=${output}.sorted.bam OUTPUT=${output}.sorted.dp.bam METRICS_FILE=${output}.sorted.dp.metrix ASSUME_SORTED=true REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true
LOG "#3 MarkDuplicates done."
  
# 4 Stat
echo "#4"
echo "${samtools} flagstat ${output}.sorted.dp.bam > ${stat}.stat"
${samtools} flagstat ${output}.sorted.dp.bam > ${stat}.stat
LOG "#4 FlagStat done."

# 5
echo "#5"
echo "java -Xmx15g -XX:ParallelGCThreads=8 -Djava.io.tmpdir=${tmp} -jar ${gatk} -T RealignerTargetCreator -nt 7 -R ${ref} -I ${output}.sorted.dp.bam -o ${output}.intervals -known ${mills} -known ${KG}"
LOG "#5 RealignerTargetCreator start."
java -Xmx15g -XX:ParallelGCThreads=8 -Djava.io.tmpdir=${tmp} -jar ${gatk} -T RealignerTargetCreator -nt 7 -R ${ref} -I ${output}.sorted.dp.bam -o ${output}.intervals -known ${mills} -known ${KG}
LOG "#5 RealignerTargetCreator done."

# 6
echo "#6"
echo "java -Xmx15g -XX:ParallelGCThreads=8 -Djava.io.tmpdir=${tmp} -jar ${gatk} -T IndelRealigner -R ${ref} -I ${output}.sorted.dp.bam -o ${output}.sorted.dp.ir.bam -targetIntervals ${output}.intervals  -known ${mills} -known ${KG}"
java -Xmx15g -XX:ParallelGCThreads=8 -Djava.io.tmpdir=${tmp} -jar ${gatk} -T IndelRealigner -R ${ref} -I ${output}.sorted.dp.bam -o ${output}.sorted.dp.ir.bam -targetIntervals ${output}.intervals  -known ${mills} -known ${KG}
LOG "#6 IndelRealigner done."

# 7
echo "#7"
echo "java -Xmx15g -XX:ParallelGCThreads=8 -Djava.io.tmpdir=${tmp} -jar ${gatk} -T BaseRecalibrator -nct 7 -I ${output}.sorted.dp.ir.bam -R ${ref} -knownSites ${dbsnp} -knownSites ${mills} -knownSites ${KG} -o ${output}.sorted.dp.ir.grp"
java -Xmx15g -XX:ParallelGCThreads=8 -Djava.io.tmpdir=${tmp} -jar ${gatk} -T BaseRecalibrator -nct 7 -I ${output}.sorted.dp.ir.bam -R ${ref} -knownSites ${dbsnp} -knownSites ${mills} -knownSites ${KG} -o ${output}.sorted.dp.ir.grp
LOG "#7 BaseRecalibrator done."

# 8
echo "#8"
echo "java -Xmx15g -XX:ParallelGCThreads=8 -Djava.io.tmpdir=${tmp} -jar ${gatk} -T PrintReads -nct 7 -R ${ref} -I ${output}.sorted.dp.ir.bam -BQSR ${output}.sorted.dp.ir.grp -o ${output}.recal.bam"
java -Xmx15g -XX:ParallelGCThreads=8 -Djava.io.tmpdir=${tmp} -jar ${gatk} -T PrintReads -nct 7 -R ${ref} -I ${output}.sorted.dp.ir.bam -BQSR ${output}.sorted.dp.ir.grp -o ${output}.recal.bam
LOG "#8 PrintReads done."

# 9 Stat
echo "#9"
echo "java -Xmx7g -XX:ParallelGCThreads=4 -Djava.io.tmpdir=${tmp} -jar ${picard_folder}CollectInsertSizeMetrics.jar MAX_RECORDS_IN_RAM=350000 VALIDATION_STRINGENCY=LENIENT INPUT=${output}.recal.bam OUTPUT=${stat}.stat.insertsizematrix HISTOGRAM_FILE=${stat}.stat.hist.pdf"
echo "java -Xmx7g -XX:ParallelGCThreads=4 -Djava.io.tmpdir=${tmp} -jar ${gatk} -T DepthOfCoverage -R ${ref} -I ${output}.recal.bam -omitBaseOutput -omitLocusTable -omitIntervals -o ${stat}.cov.out --minMappingQuality 20 --minBaseQuality 20 --summaryCoverageThreshold 8 --summaryCoverageThreshold 15 --summaryCoverageThreshold 30"
(java -Xmx7g -XX:ParallelGCThreads=4 -Djava.io.tmpdir=${tmp} -jar ${picard_folder}CollectInsertSizeMetrics.jar MAX_RECORDS_IN_RAM=350000 VALIDATION_STRINGENCY=LENIENT INPUT=${output}.recal.bam OUTPUT=${stat}.stat.insertsizematrix HISTOGRAM_FILE=${stat}.stat.hist.pdf)&
(java -Xmx7g -XX:ParallelGCThreads=4 -Djava.io.tmpdir=${tmp} -jar ${gatk} -T DepthOfCoverage -R ${ref} -I ${output}.recal.bam -omitBaseOutput -omitLocusTable -omitIntervals -o ${stat}.cov.out --minMappingQuality 20 --minBaseQuality 20 --summaryCoverageThreshold 8 --summaryCoverageThreshold 15 --summaryCoverageThreshold 30)&
wait
LOG "#9 DepthOfCoverage done."

LOG "Done sample ${sample_id}, total elasped time: $((${time_now} - ${time_start})) sec."

# cleanup
#rm $output.sam && rm $output.sorted.bam
#rm $output.sorted.dp.bam && rm $output.sorted.dp.ir.bam && rm $output.sorted.dp.ir.grp 
