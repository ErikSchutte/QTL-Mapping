##### Library ####
library(MatrixEQTL)
##### Data ####
## deliverable
load("~/Dropbox/Erik Schutte Internship 2016/Data/eQTL-data/gsTcell_genotype_frequencies_transformed.Rdata")
load("~/Dropbox/Erik Schutte Internship 2016/Data/eQTL-data/gsTcell_gene_expr_vst_88samples.Rdata")
load("~/Dropbox/Erik Schutte Internship 2016/Data/eQTL-data/gsTcell_genotypes_transformed.Rdata")
## exon count htseq
load("~/Dropbox/Erik Schutte Internship 2016/Data/eQTL-data/gsTcell_exon_expr_vst_88samples.Rdata")
load("~/Dropbox/Erik Schutte Internship 2016/Data/eQTL-data/gsTcell_genotypes_transformed_gencode.Rdata")

##### Settings ####
## Cis window
cisDist <- 5e5

## Statistical model
useModel = modelLINEAR

## Genotype file + Location file
# genotype.data <- snps.t
genotype.data <- snps.t.gencode
# genotype.loc <-read.table("~/Dropbox/Erik Schutte Internship 2016/Data/eQTL-data/snp.txt",
#                           header=T, stringsAsFactors = F)
genotype.loc <- read.table("~/Dropbox/Erik Schutte Internship 2016/Data/eQTL-data/snp.txt",
                           header=T, stringsAsFactors = F)
## Gene expression + Location file
# expression.data <- ge
expression.data <- ge.gencode
# expression.loc <- read.table("~/Dropbox/Erik Schutte Internship 2016/Data/eQTL-data/genepos.txt",
#                              header=T, stringsAsFactors = F)
expression.loc <- read.table("~/Dropbox/Erik Schutte Internship 2016/Data/eQTL-data/gencode_genes.txt",
                         stringsAsFactors = F)
## Covariates
covariates <- character() # For no covariates.
covariates <- read.table("~/Dropbox/Erik Schutte Internship 2016/Data/eQTL-data/covariates_leiden_samples.txt",
                         header=T, stringsAsFactors = F, row.names = 1)
covariates <- covariates[,mixedsort(colnames(covariates))]

## Output file names
output_cis <- "~/Dropbox/Erik Schutte Internship 2016/Results/Cis/GENCODE/new_output_cis.tsv"
output_trans <- "~/Dropbox/Erik Schutte Internship 2016/Results/Trans/GENCODE/new_output_tra.tsv" 

## Pvalue cutoffs
pvCis <- 0.5
pvTrans <- 1e-5

## Error covariance
errorCov <- numeric()

## Load genotype data
snps <- SlicedData$new();
genotype.data <- cbind(genotype.data, genotype.data, genotype.data, genotype.data)
snps <- snps$CreateFromMatrix(as.matrix(genotype.data));
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows

## Load gene expression data

gene = SlicedData$new();
gene <- gene$CreateFromMatrix(as.matrix(expression.data));
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows

## Load covariates

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name);
}

me = Matrix_eQTL_main(
  snps = snps, 
  gene = gene, 
  cvrt = cvrt,
  output_file_name     = output_trans,
  pvOutputThreshold     = pvTrans,
  useModel = useModel, 
  errorCovariance = errorCov, 
  verbose = TRUE, 
  output_file_name.cis = output_cis,
  pvOutputThreshold.cis = pvCis,
  snpspos = genotype.loc, 
  genepos = expression.loc,
  cisDist = cisDist,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE)

save(me, file="~/Dropbox/Erik Schutte Internship 2016/Data/eQTL-data/gsTcell_qtlmapping_all_timepoitns.Rdata")
