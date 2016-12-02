#!/bin/bash

SAMPLE_NAME=$1
BAM=$2
VCF=$3

module load samtools 
module load bcftools  

FASTA=/publicdata/hg19_20151104/hg19_sorted.fa
SNPBEDFILE=/frazer01/home/mdonovan/sample_identity_check/reference/1kg_AF45_55_snps.bed

samtools mpileup -D -l ${SNPBEDFILE} -uf ${FASTA} ${BAM} | bcftools call -m -Oz > ${VCF}/${SAMPLE_NAME}.vcf.gz
