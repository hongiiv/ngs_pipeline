#!/bin/bash
#$ -V
#$ -S /bin/bash

# Somatic variant call script for Tumor exsome sample
# runs on 8vCore, 16GB RAM machine

source ./pipeline.cfg

# arguments
sample_id=$1
cluster_name=$2

# log function
time_start=$(date +%s)
time_now=${time_start}
time_end=${time_start}
LOG_FILE="./logs/${cluster_name}-${sample_id}.log"
LOG() {
  time_now=$(date +%s)
  echo "[${sample_id}][$(date "+%F %T")][time: $((${time_now} - ${time_end})) sec] $*" >> ${LOG_FILE}
  time_end=${time_now}
}

# variables
bam_dir="${work_dir}results/charles/bam/"
mutect_dir="${work_dir}results/charles/mutect/"
somaticindel_dir="${work_dir}results/charles/somaticindel/"
varscan_dir="${work_dir}results/charles/varscan/"
rpkm_dir="${work_dir}results/charles/rpkm/"
tmp="/SNUH/SNUH_2/tmp/${sample_id}"

LOG "Start processing ${sample_id}..."

mkdir ${tmp}
apt-get purge -y oracle-java7-installer

#1 Matched Normal & Tumor
echo""
echo "java -Xmx16g -XX:ParallelGCThreads=8 -Djava.io.tmpdir=${tmp} -jar ${muTect} --analysis_type MuTect --reference_sequence ${ref} --dbsnp ${dbsnp_vc} --input_file:normal ${bam_dir}${sample_id}N.recal.bam --input_file:tumor ${bam_dir}${sample_id}T.recal.bam --normal_sample_name ${sample_id}N --tumor_sample_name ${sample_id}T --only_passing_calls --out ${mutect_dir}${sample_id}.muTect.call_stats.txt --vcf ${mutect_dir}${sample_id}.muTect.vcf"
java -Xmx16g -XX:ParallelGCThreads=8 -Djava.io.tmpdir=${tmp} -jar ${muTect} --analysis_type MuTect --reference_sequence ${ref} --dbsnp ${dbsnp_vc} --input_file:normal ${bam_dir}${sample_id}N.recal.bam --input_file:tumor ${bam_dir}${sample_id}T.recal.bam --normal_sample_name ${sample_id}N --tumor_sample_name ${sample_id}T --only_passing_calls --out ${mutect_dir}${sample_id}.muTect.call_stats.txt --vcf ${mutect_dir}${sample_id}.muTect.vcf
LOG "Done MuTect Call"

#2 Matched Normal & Tumor
echo""
echo "java -Xmx16g -XX:ParallelGCThreads=8 -Djava.io.tmpdir=${tmp} -jar ${indelgenotyper} -T IndelGenotyperV2 -R ${ref} -ws 300 -somatic -minCoverage 8 -minNormalCoverage 8 -I:normal ${bam_dir}${sample_id}N.recal.bam -I:tumor ${bam_dir}${sample_id}T.recal.bam -o ${somaticindel_dir}${sample_id}.somaticIndel.vcf"
java -Xmx16g -XX:ParallelGCThreads=8 -Djava.io.tmpdir=${tmp} -jar ${indelgenotyper} -T IndelGenotyperV2 -R ${ref} -ws 300 -somatic -minCoverage 8 -minNormalCoverage 8 -I:normal ${bam_dir}${sample_id}N.recal.bam -I:tumor ${bam_dir}${sample_id}T.recal.bam -o ${somaticindel_dir}${sample_id}.somaticIndel.vcf
LOG "Done somaticIndel"

#3 Matched Normal & Tumor
echo ""
echo "cp ${bam_dir}${sample_id}N.recal.bai ${bam_dir}${sample_id}N.recal.bam.bai"
echo "cp ${bam_dir}${sample_id}T.recal.bai ${bam_dir}${sample_id}T.recal.bam.bai"
echo "$samtools mpileup -Q 20 -q 20 -l ${intervals} -f ${ref} ${bam_dir}${sample_id}N.recal.bam > ${varscan_dir}${sample_id}.normal.out"
echo "$samtools mpileup -Q 20 -q 20 -l ${intervals} -f ${ref} ${bam_dir}${sample_id}T.recal.bam > ${varscan_dir}${sample_id}.tumor.out"
echo "java -Xmx16g -XX:ParallelGCThreads=8 -Djava.io.tmpdir=${tmp} -jar $varscan somatic ${varscan_dir}${sample_id}.normal.out ${varscan_dir}${sample_id}.tumor.out ${varscan_dir}${sample_id}.basename" --output-vcf 1 
cp ${bam_dir}${sample_id}N.recal.bai ${bam_dir}${sample_id}N.recal.bam.bai
cp ${bam_dir}${sample_id}T.recal.bai ${bam_dir}${sample_id}T.recal.bam.bai
$samtools mpileup -Q 20 -q 20 -l ${intervals} -f ${ref} ${bam_dir}${sample_id}N.recal.bam > ${varscan_dir}${sample_id}.normal.out
$samtools mpileup -Q 20 -q 20 -l ${intervals} -f ${ref} ${bam_dir}${sample_id}T.recal.bam > ${varscan_dir}${sample_id}.tumor.out
java -Xmx16g -XX:ParallelGCThreads=8 -Djava.io.tmpdir=${tmp} -jar $varscan somatic ${varscan_dir}${sample_id}.normal.out ${varscan_dir}${sample_id}.tumor.out ${varscan_dir}${sample_id}.basename --output-vcf 1
LOG "Done varscan"

# 4 For each sample
echo ""
echo "$python ${conifer} rpkm --probes ${intervals} --input ${bam_dir}${sample_id}N.recal.bam --output ${rpkm_dir}${sample_id}N.rpkm.txt"
echo "$python ${conifer} rpkm --probes ${intervals} --input ${bam_dir}${sample_id}T.recal.bam --output ${rpkm_dir}${sample_id}T.rpkm.txt"
$python ${conifer} rpkm --probes ${intervals} --input ${bam_dir}${sample_id}N.recal.bam --output ${rpkm_dir}${sample_id}N.rpkm.txt
$python ${conifer} rpkm --probes ${intervals} --input ${bam_dir}${sample_id}T.recal.bam --output ${rpkm_dir}${sample_id}T.rpkm.txt
LOG "Done rpkm"

# Upload

# Clean up
#rm ${varscan_dir}${sample_id}.normal.out ${varscan_dir}${sample_id}.tumor.out

LOG "Done sample ${sample_id}! Total Elapsed Time: $((${time_now} - ${time_start})) seconds."
