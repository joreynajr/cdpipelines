#!/bin/bash

sortedIN=$1
plinkOUT=$2
DATA_ID=$3

module load plink/1.90b3x
module load vcftools 

reference="/frazer01/home/mdonovan/sample_identity_check/reference/sorted_master_plinkID"
for v in ${sortedIN}/*.vcf;
	do
		sample=${v##*/}
		sample=${sample%.vcf}
		sample=${sample##*_} #data id
		out=$plinkOUT #/$sample
		if [ ! -d "$out" ]; then
			mkdir $out
		fi
		vcftools --vcf $v --minDP 10 --recode --out $out/sorted_sampleDP
		awk -v OFS="\t" '{if($0!~/^#/) { $3=$1":"$2} print $0}' $out/sorted_sampleDP.recode.vcf > $out/sorted_sample_plinkID.vcf
		plink --vcf $out/sorted_sample_plinkID.vcf --make-bed --out  $out/sorted_sample_plinkID 
		plink --bfile $out/sorted_sample_plinkID --geno 0 --make-bed --out  $out/sorted_sample_plinkID
		plink --bfile $out/sorted_sample_plinkID --bmerge $reference --make-bed --out $out/merged --geno 0
		plink --vcf $out/sorted_sample_plinkID.vcf --exclude $out/merged-merge.missnp --make-bed --out  $out/sorted_sample_plinkID
		plink --bfile $out/sorted_sample_plinkID --geno 0 --make-bed --out  $out/sorted_sample_plinkID 
		plink --bfile $out/sorted_sample_plinkID --bmerge $reference --make-bed --out $out/merged --geno 0
		plink --bfile $out/merged --genome full --out $out/merged
		egrep "FID|"$DATA_ID $out/merged.genome > $out/merged.genome.sample
done
