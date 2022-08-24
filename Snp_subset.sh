#!/bin/bash -l

vcftools --vcf WG_SNPs.vcf --positions Subset2.list --recode --recode-INFO-all --out Filtered_SNPs.vcf
