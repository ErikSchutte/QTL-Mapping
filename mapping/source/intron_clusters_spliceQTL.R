##### LIBS ###
## Matrix eQTL
# mapping splice-QTLs with differntial intron excision data as phenotype matrix for Matrix EQTL.
library(MatrixEQTL)

## Genotype data
# Load genotype data.
snps <- read.table("~/Dropbox/Erik Schutte Internship 2016/Data/filtered-data/snps.t.tsv",
                   header=T, sep=" ")

snps.loc <- read.table("~/Dropbox/Erik Schutte Internship 2016/Data/eQTL-data/snp.txt",
                       header = TRUE, stringsAsFactors = FALSE)

## Expression data
# Load intron data, scaled and normalized.
introns <- read.table("~/Dropbox/Erik Schutte Internship 2016/Data/Leafcutter/intron_counts_scaled_qnorm.tsv",
                      header=T, sep="\t")

# Filter out known faulty samples .
introns <- introns[,-grep("batch2_TCC.02.*",colnames(introns))] # -> e.g batch 2 tcc 02 only has 3 timepoints measured.
introns <- introns[,-grep("batch[0-9]+_TCC[0-9]+.*",colnames(introns))] # -> norway samples have a different format, they should be removed.
introns <- introns[,-grep("batch4_TCC.17.*", colnames(introns))] # -> batch 17 only has 3 itmepoints measured.

# Rename samples names to match samples from snp data so we can compare which are there and which arent. still 2 samples that have to be left out.
names <- sapply(colnames(introns), gsub, pattern="\\.", replacement="_", USE.NAMES = F) 
names <- sapply(names, gsub, pattern="batch[0-9]_", replacement="", USE.NAMES = F)
names <- gsub(names, pattern="_t[0-9]+", replacement="")

# Check which samples are in the snp data. (84)
introns.filtered <- introns[,which(names %in% colnames(snps))]
write.table(x = introns.filtered, file = "~/Dropbox/Erik Schutte Internship 2016/Data/Leafcutter/scaled_norm_intron_counts_filtered.csv",
            sep="\t")
colnames(introns.filtered) <- names[which(names %in% colnames(snps))] # change sample names to match snp sample names.

# Load introns.bed file.
introns.bed <- read.table("~/Dropbox/Erik Schutte Internship 2016/Data/Leafcutter/intron.bed", 
                          stringsAsFactors = F, quote = "")

#introns.test <- cbind(as.numeric(introns.bed[,1]), as.numeric(introns.bed[,2]), as.numeric(introns.bed[,3]), as.numeric(introns.bed[,4]))
##### MAP ###
## eQTL mapping
# Matrix eQTL mapping -> intron start:stop positions within 100kb of snp -> cis region.
useModel = modelLINEAR

# Only associations significant at this level will be saved.
#pvOutputThreshold_tra = 1e-5
pvOutputThreshold_cis = 0.05

# Distance for local gene-SNP pairs
cisDist = 1e5

# Set covariates
covariates_file_name = "~/Dropbox/Erik Schutte Internship 2016/Data/eQTL-data/covariates_leiden_samples.txt"

# Error covariance matrix.
errorCovariance = numeric();

# Empty variable for cis- trans- qtls.
transqtls <- NULL
cisqtls <- NULL

# Multiply genotypes from 22 to 88 samples.
snps <- cbind(snps, snps) # 44
snps <- cbind(snps, snps) # 88

# Load genotype data.
snps.sd = SlicedData$new()
snps.sd$CreateFromMatrix(as.matrix(snps))
snps.sd$fileDelimiter = "\t";      # the TAB character
snps.sd$fileOmitCharacters = "NA"; # denote missing values;
snps.sd$fileSkipRows = 1;          # one row of column labels
snps.sd$fileSkipColumns = 1;       # one column of row labels
snps.sd$fileSliceSize = 2000;      # read file in slices of 2,000 rows

# Set pattern names.
pattern_name.cis <- "~/Dropbox/Erik Schutte Internship 2016/Data/Leafcutter/intron_splice_sqtl_mapping.csv"
#pattern_name.trans <- "~/Dropbox/Erik Schutte Internship 2016/Data/Leafcutter/lncrna_interactive_eqtl_trans_all_timepoints.csv"
# 
# #transqtls <- format(transqtls,scientific=TRUE)
# #cisqtls <- format(cisqtls,scientific=TRUE)
# 
# cat("Writing to file...\n\n")
# #write.csv( transqtls,file=paste(file.path(mainDir, subDir, "Trans", fsep = ""), "/lncrna_eqtls_trans.csv", sep="") )
# #write.csv( cisqtls,file=paste(file.path(mainDir, subDir, "Cis", fsep = ""), "/lncrna_eqtls_cis.csv", sep="") )
# write.csv( MAFGenotype, file=paste(file.path(mainDir, subDir, fsep = ""), "/MAFGenotype.csv", sep="")  )
# 
# cat("Wrote file(s) to: \n",paste(file.path(mainDir, subDir, "Trans", fsep=""), sep=""),
#     paste(file.path(mainDir, subDir, "Cis", fsep="") ) )
# 
# cat("\n\nFINISHED\n\n")

# Create output files.
file.create(pattern_name.cis)
#file.create(pattern_name.trans)

# Output file name and location.
output_file_name_cis = pattern_name.cis
#output_file_name_tra = pattern_name.trans

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
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
cvrt$LoadFile(covariates_file_name);

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

# Unlink output_file_name (keep it commented to keep eQTL mappings),
#unlink(output_file_name_cis)
#unlink(output_file_name_tra)

# Results.
cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n\n\n')

# Verbose - Show all detected eQTLs.
if ( VERBOSE == TRUE ) {
  cat('Detected eQTLs:', '\n')
  show(me$all$neqtls)
}