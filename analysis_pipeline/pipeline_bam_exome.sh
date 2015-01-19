#!/bin/bash
#$ -V
#$ -S /bin/bash

source ./pipeline.cfg

# arguments
file_name=$1

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
fastq_r1="/BIO/exome/T2D_family/${file_name}_1.fastq.gz"
fastq_r2="/BIO/exome/T2D_family/${file_name}_2.fastq.gz"
output="${work_dir}/bam/${file_name}"
fastqc="${work_dir}/fastqc/"

LOG "Start processing sample ${file_name}..."

# 0 Stat
(/dist_vol/app/FastQC/fastqc -t 4 -o ${fastqc} -f fastq ${fastq_r1})&
(/dist_vol/app/FastQC/fastqc -t 4 -o ${fastqc} -f fastq ${fastq_r2})&
wait
LOG "#0 FastQC done."

# 1
echo "#1"
echo "${bwa} mem -M -t 8 ${ref} ${fastq_r1} ${fastq_r2} > ${output}.sam"
${bwa} mem -M -t 8 ${ref} ${fastq_r1} ${fastq_r2} > ${output}.sam
LOG "#1 BWA done."

# 2
echo "#2"
echo "java -Xmx15g -XX:ParallelGCThreads=8 -Djava.io.tmpdir=${tmp} -jar ${picard_folder}AddOrReplaceReadGroups.jar INPUT=${output}.sam OUTPUT=${output}.sorted.bam SORT_ORDER=coordinate RGID=${file_name} RGLB=${file_name} RGPL=illumina RGPU=SureSelectAllExon RGSM=${file_name} VALIDATION_STRINGENCY=LENIENT"
java -Xmx15g -XX:ParallelGCThreads=8 -Djava.io.tmpdir=${tmp} -jar ${picard_folder}AddOrReplaceReadGroups.jar INPUT=${output}.sam OUTPUT=${output}.sorted.bam SORT_ORDER=coordinate RGID=${file_name} RGLB=${file_name} RGPL=illumina RGPU=SureSelectAllExon RGSM=${file_name} VALIDATION_STRINGENCY=LENIENT
LOG "#2 AddOrReplaceReadGroups done."
#check_line_size_bam ${fastq_r1} ${output}.sorted.bam

# 3
echo "#3"
echo "java -Xmx15g -XX:ParallelGCThreads=8 -Djava.io.tmpdir=${tmp} -jar ${picard_folder}MarkDuplicates.jar INPUT=${output}.sorted.bam OUTPUT=${output}.sorted.dp.bam METRICS_FILE=${output}.sorted.dp.metrix ASSUME_SORTED=true REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true"
java -Xmx15g -XX:ParallelGCThreads=8 -Djava.io.tmpdir=${tmp} -jar ${picard_folder}MarkDuplicates.jar INPUT=${output}.sorted.bam OUTPUT=${output}.sorted.dp.bam METRICS_FILE=${output}.sorted.dp.metrix ASSUME_SORTED=true REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT
LOG "#3 MarkDuplicates done."
  
# 4 
echo "#4"
echo "java -Xmx15g -XX:ParallelGCThreads=8 -Djava.io.tmpdir=${tmp} -jar ${picard_folder}FixMateInformation.jar INPUT=${output}.sorted.dp.bam OUTPUT=${output}.sorted.dp.fx.bam SO=coordinate VALIDATION_STRINGENCY=LENIENT"
java -Xmx15g -XX:ParallelGCThreads=8 -Djava.io.tmpdir=${tmp} -jar ${picard_folder}FixMateInformation.jar INPUT=${output}.sorted.dp.bam OUTPUT=${output}.sorted.dp.fx.bam SO=coordinate VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true
LOG "#4 FixMateInformation done."

# 5
echo "#5"
echo "java -Xmx15g -XX:ParallelGCThreads=8 -Djava.io.tmpdir=${tmp} -jar ${gatk} -T RealignerTargetCreator -nt 8 -R ${ref} -I ${output}.sorted.dp.fx.bam -o ${output}.intervals -known ${KG}"
java -Xmx15g -XX:ParallelGCThreads=8 -Djava.io.tmpdir=${tmp} -jar ${gatk} -T RealignerTargetCreator -nt 8 -R ${ref} -I ${output}.sorted.dp.fx.bam -o ${output}.intervals -known ${KG}
LOG "#5 RealignerTargetCreator done."

# 6
echo "#6"
echo "java -Xmx15g -XX:ParallelGCThreads=8 -Djava.io.tmpdir=${tmp} -jar ${gatk} -T IndelRealigner -R ${ref} -I ${output}.sorted.dp.fx.bam -o ${output}.sorted.dp.ir.bam -targetIntervals ${output}.intervals -known ${KG}"
java -Xmx15g -XX:ParallelGCThreads=8 -Djava.io.tmpdir=${tmp} -jar ${gatk} -T IndelRealigner -R ${ref} -I ${output}.sorted.dp.fx.bam -o ${output}.sorted.dp.ir.bam -targetIntervals ${output}.intervals -known ${KG}
LOG "#6 IndelRealigner done."

# 7
echo "#7"
echo "java -Xmx15g -XX:ParallelGCThreads=8 -Djava.io.tmpdir=${tmp} -jar ${gatk} -T BaseRecalibrator -nct 8 -I ${output}.sorted.dp.ir.bam -R ${ref} -knownSites ${dbsnp} -knownSites ${KG} -o ${output}.sorted.dp.ir.grp"
java -Xmx15g -XX:ParallelGCThreads=8 -Djava.io.tmpdir=${tmp} -jar ${gatk} -T BaseRecalibrator -nct 8 -I ${output}.sorted.dp.ir.bam -R ${ref} -knownSites ${dbsnp} -knownSites ${KG} -o ${output}.sorted.dp.ir.grp
LOG "#7 BaseRecalibrator done."

# 8
echo "#8"
echo "java -Xmx15g -XX:ParallelGCThreads=8 -Djava.io.tmpdir=${tmp} -jar ${gatk} -T PrintReads -nct 8 -R ${ref} -I ${output}.sorted.dp.ir.bam -BQSR ${output}.sorted.dp.ir.grp -o ${output}.recal.bam"
java -Xmx15g -XX:ParallelGCThreads=8 -Djava.io.tmpdir=${tmp} -jar ${gatk} -T PrintReads -nct 8 -R ${ref} -I ${output}.sorted.dp.ir.bam -BQSR ${output}.sorted.dp.ir.grp -o ${output}.recal.bam
LOG "#8 PrintReads done."

LOG "Done sample ${file_name}, total elasped time: $((${time_now} - ${time_start})) sec."
