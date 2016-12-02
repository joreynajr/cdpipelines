#!/bin/bash

#$ -N sampleMatchingRNA
#$ -pe smp 4      ### specify number of cores requested
#$ -l week       ### specify queue ([short, week, long, opt], default is all)
#$ -e /frazer01/home/mdonovan/twin_variation/rna_analysis/variants/sh/run.out   ### redirect stderr to this file
#$ -o /frazer01/home/mdonovan/twin_variation/rna_analysis/variants/sh/run.err   ### redirect stdout to this file

# Module load cardips should override this.
module load cardips  

FASTA=/publicdata/hg19_20151104/hg19_sorted.fa
SNPBEDFILE=/frazer01/home/mdonovan/sample_identity_check/reference/1kg_AF45_55_snps.bed

if [ $# -eq 0 ]
  then
    echo "ERROR in RNA Sample Match Pipeline. Please supply path to data ID."
fi

DATA_ID_PATH=$1 #user provides path to data id, e.g. /projects/CARDIPS/pipeline/RNAseq/sample/data_id
DATA_ID=${DATA_ID_PATH##*/}
DATA_ID_BAM=$DATA_ID_PATH/alignment/$DATA_ID"_sorted_mdup.bam"

vcf_folder=$DATA_ID_PATH/qc/$DATA_ID"_sample_swap"
PLINK=${vcf_folder}"/plink"

if [ ! -d "$vcf_folder" ]; then
	mkdir $vcf_folder
fi

if [ ! -d "$PLINK" ]; then
	mkdir $PLINK
fi

echo "calling variants on: " $DATA_ID
variant_file=${vcf_folder}/${DATA_ID}"_1kg_variants.vcf.gz"
source callVariants.sh $DATA_ID $DATA_ID_BAM $vcf_folder
source modifyVCFs.sh $variant_file $vcf_folder
source plink.sh $vcf_folder $PLINK $DATA_ID
