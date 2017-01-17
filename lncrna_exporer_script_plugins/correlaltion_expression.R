##########################################################
#Author: Marije van der Geest & Erik Schutte             #
#                                                        #
#Generates line plot for gluten specific T-cell          #
#expression per gene under different conditions          #
#(timpoints 0, 10, 30, 180)                              #
##########################################################

### Set-up
## Libraries.
# Libraries used to create the plot.
library(gplots)
library(corrplot)
library(RColorBrewer)

## System login.
# Login to the Molgenis R API.
source("http://localhost:8080/molgenis.R")
molgenis.login("admin", "admin")

### Functions
## Correlation matrix
# Calculate the p-value for the correlation matrix.
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

### Data
## Get selected data by user.
# Save the genes in a variable and pre-prepare them for the query.
genes <- "ENSG00000019144,ENSG00000118094,ENSG00000110367,ENSG00000160683,ENSG00000186174,ENSG00000110375,ENSG00000176302,ENSG00000186166,ENSG00000118181,ENSG00000196655,ENSG00000149428,ENSG00000256269,ENSG00000188486,ENSG00000172269,ENSG00000172375,ENSG00000255176,ENSG00000255239,ENSG00000207462,ENSG00000255422,ENSG00000254478,ENSG00000254621,ENSG00000201535,ENSG00000245869,ENSG00000264211,ENSG00000239726,ENSG00000264523,ENSG00000222529,ENSG00000242712,ENSG00000255121,ENSG00000240970,ENSG00000254428,ENSG00000254909,ENSG00000266398,ENSG00000137700,ENSG00000255114,ENSG00000271751,ENSG00000272186,ENSG00000160695,"
genes <- unlist(strsplit(genes, '[,]'))

## Prepared Statement.
# Create query for data filtration.
qs <- ""
for (i in 1:length(genes)){
  qs[i] <- paste("EnsemblGeneID==", genes[i], sep="")
}
q=paste(qs, collapse=",")

## Get data
# Load the data from the Molgenis Repository using the Prepared statement.
LLDeep <- molgenis.get("lncrna_data_Tcell123samples", num=100000, q=q)

# Set the first column from the table as row names.
rownames(LLDeep) <- LLDeep[,1]

# Remove the first column.
LLDeep <- LLDeep[,-1]

# New prepared statement for gene name conversion.
qs <- ""
for (i in 1:length(genes)){
  qs[i] <- paste("EnsemblGeneID==", genes[i], sep="")
}
q=paste(qs, collapse=",")

### Data Pre-processing
## Create time indices for the data.
# Save the index on the dataset where the data changes per time point. So all expression over time 0, 10 ,30 and 180.
timepoints <- c("t0", "t10", "t30", "t180")
index <- lapply(timepoints, function(x) {which(lapply(colnames(LLDeep), function(y) grep(y, pattern=x) ) == 1 ) } )

## Get the Ensemble Gene IDs and the associated gene names.
# Save the Information about the Ensemble Genes, their associated gene nemes etc.
conversion.table <- molgenis.get("lncrna_GeneInfo", q=q)

# Change the AUTO ID to Ensemble Gene IDs.
rownames(conversion.table) <- conversion.table[,1]

### Data Plot
## Create correlation matrix.
# Create image
#png("${outputFile}")
if ( length(genes) < 10 ) {
  pixels = 500;
} else if ( length(genes ) < 20 ) {
  pixels = 1000;
} else if ( length (genes) < 30 ) {
  pixels = 1500;
} else {
  pixels = 1600;
}
par(mfrow=c(2,2))

# For each timepoint, calculate the correlation.
for ( t in 1:length(timepoints) ) {
  
  #Perform correlation (Spearman)
  correlation.matrix <- cor(t(LLDeep[,index[[t]]]), method="spearman", use="pairwise.complete.obs")
  
  # Ommit NA's completly to prevent division zero errors.
  correlation.matrix <- correlation.matrix[rowSums(!is.na(correlation.matrix))!=0, colSums(!is.na(correlation.matrix))!=0]
  
  # Change the row and column names to Gene names instead of Ensemble Gene IDs
  colnames(correlation.matrix) <- rownames(correlation.matrix) <- conversion.table[which(rownames(conversion.table) %in% rownames(correlation.matrix) == T),2]
  
  # Plot the correlation plot.
  corrplot(correlation.matrix,type="upper", method="color", col=colorRampPalette(brewer.pal(9, name="RdBu"))(100),
           addCoef.col = "black", tl.col="black", tl.cex=0.6, tl.srt = 45, diag = TRUE, addgrid.col = NA,
           na.label="NA",number.cex=0.45)
  mtext(text=paste(as.character(timepoints[t])),side=1, line=4,c)
}
