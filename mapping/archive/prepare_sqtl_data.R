##### Data Prep 
#### Load Data. 
### Genotype data. ####
## Load genotype data.
snps <- read.table("~/Dropbox/Erik Schutte Internship 2016/Data/filtered-data/snps.t.tsv",
                   header=T, sep=" ")
snps.geno <- read.table("~/Dropbox/Erik Schutte Internship 2016/Data/filtered-data/snps.genotype.t.tsv",header=T)
alleles <- data.frame(alleles=snps.geno[,22],snp=rownames(snps.geno))

### Expression data. ####
## Load intron data, scaled and normalized.
introns <- read.table("~/Dropbox/Erik Schutte Internship 2016/Data/Leafcutter/intron_counts_scaled_qnorm.tsv",
                      header=T, sep="\t")

### Covariates ####
## Load covariates, age, gender and time.
covariates = read.table("~/Dropbox/Erik Schutte Internship 2016/Data/eQTL-data/covariates_leiden_samples.txt",
                        header=T, row.names = 1)

#### Prepare Data.
### Genotypes. ####
## Order genotype sample names as template for all sample names.
snps <- t(snps)[order(rownames(t(snps))),]
snps.geno <- t(snps.geno)[order(rownames(t(snps.geno))),]

### Expression. ####
## Match Expression samples as genotypes samples.
names.org <- sapply(colnames(introns), gsub, pattern="\\.", replacement="_", USE.NAMES = F) 
names.t <- sapply(names.org, gsub, pattern="batch[0-9]_", replacement="", USE.NAMES = F)
names <- sapply(names.t, gsub, pattern="_t[0-9]+", replacement="", USE.NAMES = F)

## Check which samples are in the snp data. (84)
introns.filtered <- introns[,which(names %in% rownames(snps))]
colnames(introns.filtered) <- names.t[which(names %in% rownames(snps))] # change sample names to match snp sample names.

## Order samples as SNP samples.
name.order <- order(colnames(introns.filtered))
introns.filtered <- introns.filtered[,colnames(introns.filtered)[name.order]]

## Save time index.
times <- list(t0=grep("_t0",colnames(introns.filtered)),
              t10=grep("_t10",colnames(introns.filtered)),
              t30=grep("_t30",colnames(introns.filtered)),
              t180=grep("_t180",colnames(introns.filtered)))

## Remove the time flag from the sample name.
colnames(introns.filtered) <- gsub(colnames(introns.filtered), pattern="_t[0-9]+", replacement="")

### Covariates ####
names <- gsub(colnames(covariates), pattern="\\.[0-9]", replacement = "")
covariates <- covariates[,order(names)]
colnames(covariates) <- gsub(colnames(covariates), pattern="\\.[0-9]", replacement="")
save(alleles, times, snps, snps.geno, introns.filtered, covariates, file="~/Dropbox/Erik Schutte Internship 2016/Data/Leafcutter/ordered_sqtl_data.Rdata")

