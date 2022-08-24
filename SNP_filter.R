library(vcfR)

vcf <- read.vcfR("WG_SNPs.vcf", verbose = FALSE)
localchromR <- create.chromR(vcf)

#Access specific fields of interest from the VCF file
AD <- extract.gt(vcf, element="AD", as.numeric = F, IDtoRowNames = F, extract = TRUE)
C_name<- getCHROM(vcf)
C_POS<- getPOS(vcf)
AF <- extract.info(vcf, element="AF", as.numeric = T, mask = FALSE)

#Access AD values in the AD field to compute number of reads 
AD1 <- masplit(AD, delim = ",",  sort = 0, record = 1)
AD2 <- masplit(AD, delim = ",",  sort = 0,record = 2)
NR <- ad1 + ad2

#Find positions with less than 5 reads and replace them with zero
NR[NR < 5] <- 0
(out1<- which(NR == 0))
ad1<-replace(ad1, out1, 0)
ad2<-replace(ad2, out1, 0)

#Compute allele fractions based on the Alt and Ref values in the AD fields
freq <- ad2/(ad1+ad2)
Af<-round(freq, digits = 3)

#Mark NaNs as missing data
Af[is.na(Af)] <-"."

#Assign sites with Allele fractions of more than 0.1 and less than 0.9 as missing data
Af[Af > 0.1 & Af < 0.9] <- "."

#Assign sites with Allele fractions of at least 0.9 as Alt alleles (1)
Af[Af >=  0.9] <- 1

#Assign sites with Allele fractions of 0.1 and less as REF alleles (0)
Af[Af <= 0.1 & Af != "."] <- 0

#Combine Chromosome Names, Chromosome positions, and the computed Allele fractions into a dataframe
File <- cbind(C_name, C_POS, Af)

#Filter the dataframe by summing the number of missing data. Keep only rows/sites with not more than 8 out of the total 24 samples with missing data. i.e only keep sites with at least 67% data
AF <- File[rowSums(File == ".") <= 8, ]

#Generate a subset file (To be used to subset the original VCF file to retain the filtered SNPs only in a new VCF file for downstream analyses) 
Final_file<- data.frame(AF)
subset_file <- Final_file[,1:2]
write.table(subset_file, file = "vcf_subset_file.txt", quote=FALSE, col.names = FALSE, row.names = F)
