#!/bin/bash
#$ -V
#$ -S /bin/bash

# Variant calling script for WGS sample
# runs on 8vCore, 64GB RAM machine

source ./pipeline.cfg

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
for ((i=1;i<38;i++)); do
  inputs="${inputs} -V ${work_dir}results/kr37/germline/KR$i.hc.vcf"
done

mkdir ${tmp}

LOG "Start processing KR37wgs..."

## Genotype GVCFs (Use All hc.vcf files)
echo "java -Xmx60g -XX:ParallelGCThreads=8 -Djava.io.tmpdir=${tmp} -jar ${gatk} -T GenotypeGVCFs -nt 7 -R ${ref} ${inputs} -o ${out_dir}KR37wgs.vcf"
java -Xmx60g -XX:ParallelGCThreads=8 -Djava.io.tmpdir=${tmp} -jar ${gatk} -T GenotypeGVCFs -nt 7 -R ${ref} ${inputs} -o ${out_dir}KR37wgs.vcf
echo "java -Xmx60g -XX:ParallelGCThreads=8 -Djava.io.tmpdir=${tmp} -jar ${gatk} -T SelectVariants -R ${ref} --variant ${out_dir}KR37wgs.vcf --restrictAllelesTo BIALLELIC -o ${out_dir}KR37wgs.bi.vcf"
java -Xmx60g -XX:ParallelGCThreads=8 -Djava.io.tmpdir=${tmp} -jar ${gatk} -T SelectVariants -R ${ref} --variant ${out_dir}KR37wgs.vcf --restrictAllelesTo BIALLELIC -o ${out_dir}KR37wgs.bi.vcf
echo "java -Xmx60g -XX:ParallelGCThreads=8 -Djava.io.tmpdir=${tmp} -jar ${gatk} -T SelectVariants -R ${ref} -V ${out_dir}KR37wgs.bi.vcf -o ${out_dir}KR37wgs.snp.vcf -selectType SNP"
java -Xmx60g -XX:ParallelGCThreads=8 -Djava.io.tmpdir=${tmp} -jar ${gatk} -T SelectVariants -R ${ref} -V ${out_dir}KR37wgs.bi.vcf -o ${out_dir}KR37wgs.snp.vcf -selectType SNP
echo "java -Xmx60g -XX:ParallelGCThreads=8 -Djava.io.tmpdir=${tmp} -jar ${gatk} -T SelectVariants -R ${ref} -V ${out_dir}KR37wgs.bi.vcf -o ${out_dir}KR37wgs.indel.vcf -selectType INDEL"
java -Xmx60g -XX:ParallelGCThreads=8 -Djava.io.tmpdir=${tmp} -jar ${gatk} -T SelectVariants -R ${ref} -V ${out_dir}KR37wgs.bi.vcf -o ${out_dir}KR37wgs.indel.vcf -selectType INDEL

## Variant Filtration (SNP)
echo "java -Xmx60g -XX:ParallelGCThreads=8 -Djava.io.tmpdir=${tmp} -jar ${gatk} -T VariantFiltration -R ${ref} -V ${out_dir}KR37wgs.snp.vcf -o ${out_dir}KR37wgs.snp.filtered.vcf -window 35 -cluster 3 --filterName FS -filter 'FS > 30.0' --filterName LowQD -filter 'QD < 2.0'"
echo "java -Xmx60g -XX:ParallelGCThreads=8 -Djava.io.tmpdir=${tmp} -jar ${gatk} -T VariantFiltration -R ${ref} -V ${out_dir}KR37wgs.indel.vcf -o ${out_dir}KR37wgs.indel.filtered.vcf -window 35 -cluster 3 -filterName FS -filter 'FS > 30.0' --filterName QD -filter 'QD < 2.0'"
java -Xmx60g -XX:ParallelGCThreads=8 -Djava.io.tmpdir=${tmp} -jar ${gatk} -T VariantFiltration -R ${ref} -V ${out_dir}KR37wgs.snp.vcf -o ${out_dir}KR37wgs.snp.filtered.vcf -window 35 -cluster 3 --filterName FS -filter "FS > 30.0" --filterName LowQD -filter "QD < 2.0"
java -Xmx60g -XX:ParallelGCThreads=8 -Djava.io.tmpdir=${tmp} -jar ${gatk} -T VariantFiltration -R ${ref} -V ${out_dir}KR37wgs.indel.vcf -o ${out_dir}KR37wgs.indel.filtered.vcf -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" --filterName QD -filter "QD < 2.0"

## Select Unfiltered Variants only (SNP)
echo "java -Xmx60g -XX:ParallelGCThreads=8 -Djava.io.tmpdir=${tmp} -jar ${gatk} -T SelectVariants -R ${ref} -V ${out_dir}KR37wgs.snp.filtered.vcf -o ${out_dir}KR37wgs.snp.filtered_only.vcf -env -ef"
echo "java -Xmx60g -XX:ParallelGCThreads=8 -Djava.io.tmpdir=${tmp} -jar ${gatk} -T SelectVariants -R ${ref} -V ${out_dir}KR37wgs.indel.filtered.vcf -o ${out_dir}KR37wgs.indel.filtered_only.vcf -env -ef"
java -Xmx60g -XX:ParallelGCThreads=8 -Djava.io.tmpdir=${tmp} -jar ${gatk} -T SelectVariants -R ${ref} -V ${out_dir}KR37wgs.snp.filtered.vcf -o ${out_dir}KR37wgs.snp.filtered_only.vcf -env -ef
java -Xmx60g -XX:ParallelGCThreads=8 -Djava.io.tmpdir=${tmp} -jar ${gatk} -T SelectVariants -R ${ref} -V ${out_dir}KR37wgs.indel.filtered.vcf -o ${out_dir}KR37wgs.indel.filtered_only.vcf -env -ef

## Merge
echo "java -Xmx60g -XX:ParallelGCThreads=8 -Djava.io.tmpdir=${tmp} -jar ${gatk} -T CombineVariants -R ${ref} --variant ${out_dir}KR37wgs.snp.filtered_only.vcf --variant ${out_dir}KR37wgs.indel.filtered_only.vcf -o ${out_dir}KR37wgs.combined.filtered_only.vcf"
java -Xmx60g -XX:ParallelGCThreads=8 -Djava.io.tmpdir=${tmp} -jar ${gatk} -T CombineVariants -R ${ref} --variant ${out_dir}KR37wgs.snp.filtered_only.vcf --variant ${out_dir}KR37wgs.indel.filtered_only.vcf -o ${out_dir}KR37wgs.combined.filtered_only.vcf

## Annovar
echo "perl ${annovar}convert2annovar.pl ${out_dir}KR37wgs.combined.filtered_only.vcf --format vcf4old -includeinfo --outfile ${out_dir}KR37wgs.combined.filtered_only.annovar"
echo "perl ${annovar}table_annovar.pl ${out_dir}KR37wgs.combined.filtered_only.annovar ${humandb} --buildver hg19 --protocol refGene,1000g2012apr_all,1000g2012apr_asn,gwasCatalog,cosmic67wgs,snp138NonFlagged,genomicSuperDups --operation g,f,f,r,f,f,r --outfile ${out_dir}KR37wgs.combined.filtered_only --remove -otherinfo"
perl ${annovar}convert2annovar.pl ${out_dir}KR37wgs.combined.filtered_only.vcf --format vcf4old -includeinfo --outfile ${out_dir}KR37wgs.combined.filtered_only.annovar
perl ${annovar}table_annovar.pl ${out_dir}KR37wgs.combined.filtered_only.annovar ${humandb} --buildver hg19 --protocol refGene,1000g2012apr_all,1000g2012apr_asn,gwasCatalog,cosmic67wgs,snp138NonFlagged,genomicSuperDups --operation g,f,f,r,f,f,r --outfile ${out_dir}KR37wgs.combined.filtered_only --remove -otherinfo

## Clean up
#rm ${out_dir}*.hc.vcf
#(rm KR37wgs.bi.vcf) &&(rm KR37wgs.snp.vcf) && (rm KR37wgs.indel.vcf)
#(rm KR37wgs.*filtered.vcf) &&(rm KR37wgs.snp.filtered_only.vcf) && (rm KR37wgs.indel.filtered_only.vcf)
