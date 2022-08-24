#!/bin/bash -l

vcftools --vcf WG_SNPs.vcf --positions vcf_subset_file.txt --recode --recode-INFO-all --out Filtered_SNPs.vcf
