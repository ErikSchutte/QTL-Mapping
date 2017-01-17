##########################################################
#Author: Marije van der Geest & Erik Schutte             #
#                                                        #
# This data is filtered for Norway data and NA's.        #
#                                                        #
#Generates line plot for gluten specific T-cell          #
#expression per gene under different conditions          #
#(timpoints 0, 10, 30, 180)                              #
##########################################################

### Set-up
## Libraries needed.
# grDevices.
library(grDevices)
library(ggplot2)

## Molgenis connection.
# Login to molgenis as admin admin on localhost to connect with the R-Api.
source("http://localhost:8080/molgenis.R")
molgenis.login("admin", "admin")

### Functions
## Function for mean per timepoint calculation.
# Function to calculate the means per time points.
calculateMeanPerTimepoint <- function(exp.values, times) {
  
  # For each timepoint.. 
  meanExp <- lapply(times, function(x) {
    # ..get the samples with timepoint x and calculate the mean of those values.
    mean(as.numeric(exp.values[,grep(x, names(exp.values), value=T)]))
   }
  )
  maxExp <- pmax(meanExp)
  minExp <- pmin(meanExp)
}

### Pre-processing
## Prepare data.
# Set time points for filtration.
time.labels <- c(0,10,30,180)

# Calculate the mean expression per gene for each timepoint and store it in a matrix.
times <- c("__t0","__t10","__t30","__t180")

## Plot y-axis
# To dynamically increase/decrease the y-axis of the plot save the minimum and maximum values.
maxExpressed <- -Inf
minExpressed <- Inf

## Get genes.
# Get genes from user input.
genes <- "ENSG00000019144,ENSG00000118094,ENSG00000110367,"

# Split on the , and create a string from the list.
genes <- unlist(strsplit(genes, '[,]'))

## Prepared-statement
# Create a prepared query statement to get gene expression from the database.
qs <- ""

# For all genes.
for (i in 1:length(genes)){
  # Create a query string with 'EnsemblGeneID==ENSG00000204913'.
  qs[i] <- paste("EnsemblGeneID==", genes[i], sep="")
}
# Create the query.
q=paste(qs, collapse=",")

## Retrieve data.
# Get data from the database using the prepared query.
gs.tcell <- molgenis.get("lncrna_data_Tcell123samples", q=q)

## Format data.
# Set the rownames instead of the auto numbering to ENSG's.
rownames(gs.tcell) <- gs.tcell[,1]

# Remove this row from the dataset.
gs.tcell <- gs.tcell[,-1]

### Processing
## Matrix
# Set empty matrix for all mean Expressions per timepoint.
meanExp.perTimepoint.matrix <- NULL

# For all genes in the selected genes.
for (i in 1:length(genes)) {
  
  # Subste the expression values for the selected genes.
  exp.values <- gs.tcell[genes,]
  
  # Call 'calculateMeanPerTimepoint' to calculate the mean expression per timepoint per gene.
  meanExp.perTimepoint <- calculateMeanPerTimepoint(exp.values[i,], times)
  
  # Bind the means per timepoint to the previously empty created matrix.
  meanExp.perTimepoint.matrix <- rbind(meanExp.perTimepoint.matrix, unlist(meanExp.perTimepoint))
  
  # Find the highest value for the y-axis of the plot.
  if (max(meanExp.perTimepoint.matrix[i,]) > maxExpressed) {
    maxExpressed = max(meanExp.perTimepoint.matrix[i,])
  }
  
  # Find the lowest value for the y-axis of the plot.
  if (min(meanExp.perTimepoint.matrix[i,]) < minExpressed) {
    minExpressed = min(meanExp.perTimepoint.matrix[i,])
  }
  
}

# Set the column names for the empty matrix.
colnames(meanExp.perTimepoint.matrix) <- c(0,10,30,180)

# Set the row names for the empty matrix.
rownames(meanExp.perTimepoint.matrix) <- genes

## Generate Image.
# Set color spectrum.
colors <- rainbow(length(genes))

# Add extra space to right of plot area; change clipping to figure.
frame()
par(mar=c(5.1, 4.1, 4.1, 9.1))

## Image
# Create output image.
#png("${outputFile}")
#dev.copy(png,'myplot.png')

# Create empty plot.
plot(0,type="b", col="red", xlab="Time [min]", ylab="Mean expression", xaxt="n", ylim=c((minExpressed-1),(maxExpressed+1)), xlim=c(1,4),cex.main=2.2)
axis(1, at=1:length(time.labels), labels=as.vector(time.labels))

# Draw lines.
for ( i in 1:length(genes) ) {
  lines(meanExp.perTimepoint.matrix[i,], type="b", col=colors[i], lwd=2, pch=16)
}

# Create query for gene id -> gene name conversion.
qs <- ""
for (i in 1:length(genes)){
  qs[i] <- paste("EnsemblGeneID==", genes[i], sep="")
}
q=paste(qs, collapse=",")

# Convert genes.
conversion.table <- molgenis.get("lncrna_GeneInfo", q=q)
genes = as.character(conversion.table[do.call(c,lapply(genes, function(x) which(as.character(conversion.table[,1]) == x))),2])

# Generate legend.
legend("top", legend=genes, col=colors, pch=19, cex=0.8,bty='n', horiz=F,xpd=T,
       text.width=0.4)

#dev.off()

