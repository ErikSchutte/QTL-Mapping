##### LIBS ###
## Matrix eQTL
# mapping splice-QTLs with differntial intron excision data as phenotype matrix for Matrix EQTL.
library(MatrixEQTL)

load("~/Dropbox/Erik Schutte Internship 2016/Data/Leafcutter/ordered_sqtl_data.Rdata")

# Load introns.bed file.
introns.bed <- read.table("~/Dropbox/Erik Schutte Internship 2016/Data/Leafcutter/intron.bed", 
                          stringsAsFactors = F, quote = "")

## sQTL mapping
useModel = modelLINEAR

# Only associations significant at this level will be saved.
pvOutputThreshold_cis = 0.05

# Distance for local gene-SNP pairs
cisDist = 5e5

# Set covariates
covariates_file_name = covariates

# Error covariance matrix.
errorCovariance = numeric();

snps <- cbind(t(snps), t(snps), t(snps), t(snps))
snps <- snps[,names(introns.filtered)]

# Load genotype data.
snps.sd = SlicedData$new()
snps.sd$CreateFromMatrix(snps)
snps.sd$fileDelimiter = "\t";      # the TAB character
snps.sd$fileOmitCharacters = "NA"; # denote missing values;
snps.sd$fileSkipRows = 1;          # one row of column labels
snps.sd$fileSkipColumns = 1;       # one column of row labels
snps.sd$fileSliceSize = 2000;      # read file in slices of 2,000 rows

# Set pattern names.
pattern_name.cis <- "~/Dropbox/Erik Schutte Internship 2016/Data/Leafcutter/intron_splice_sqtl_mapping.csv"

# Create output files.
file.create(pattern_name.cis)

# Output file name and location.
output_file_name_cis = pattern_name.cis

# Load gene expression data .
gene = SlicedData$new();
gene$CreateFromMatrix(as.matrix(introns.filtered))

gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows

# Load covariates.
cvrt = SlicedData$new();
cvrt$CreateFromMatrix(as.matrix(covariates_file_name));
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels


## Run the analysis.
me = Matrix_eQTL_main(
  snps = snps.sd,
  gene = gene,
  cvrt = cvrt,
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  output_file_name = output_file_name_cis,
  pvOutputThreshold = pvOutputThreshold_cis,
  snpspos = snps.loc,
  genepos = as.data.frame(introns.test),
  cisDist = cisDist,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);

# Results.
cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n\n\n')

save(me, file="~/Dropbox/Erik Schutte Internship 2016/Data/Leafcutter/intron_splice_sqtl_mapping.csv")
