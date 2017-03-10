#!/bin/bash

variant_file=$1
VCF=$2
filename=${variant_file##*/}
name=${filename%.*}

for v in ${VCF}/*.gz; 
	do
		gunzip $v;
	done

for v in ${VCF}/*.vcf; 
	do
		sed -e 's!chr!!' $v > ${VCF}/"nochr_"${name}
		rm $v
	done

for v in ${VCF}/nochr_*.vcf; 
	do
		(grep ^# $v; grep -v ^# $v|sort -k1,1n -k2,2n) > ${VCF}/${name}
		rm $v
	done

