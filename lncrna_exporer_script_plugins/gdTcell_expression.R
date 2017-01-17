##########################################################
#Authors: Marije van der Geest & Erik Schutte            #
#                                                        #
#Generates heatmap for gamma/delta T-cell expression     #
#per gene under different conditio                       #
##########################################################

### Set-up
## Libraries needed
# ggplot2.
library("gplots")
library("fields")
library("ggplot2")

# RColorBrewer.
library("RColorBrewer")

# reshape2.
require("reshape2")

## Molgenis connection.
# Login to molgenis as admin admin on localhost to connect with the R-Api.
source("http://localhost:8080/molgenis.R")
molgenis.login("admin", "admin")

### Pre-processing
## Get genes.
# Get genes from user input.
genes <- "ENSG00000204913,ENSG00000172057,ENSG00000073605,ENSG00000008838,"

# Split on the , and create a string from the list.
genes <- unlist(strsplit(genes, '[,]'))

## Prepared-statement
# Create a prepared query statement to get gene expression from the database.
qs <- ""

# For all genes.
for (i in 1:length(genes)){
  # Create a query string with e.g.'EnsemblGeneID==ENSG00000204913'.
  qs[i] <- paste("EnsemblGeneID==", genes[i], sep="")
}
# Create the query.
q=paste(qs, collapse=",")

## Retrieve data.
# Get data from the database using the prepared query.
gdTcells.means <- molgenis.get("lncrna_data_meangdTcells", q=q)

## Format data.
# Set the rownames instead of the auto numbering to ENSG's.
rownames(gdTcells.means) <- gdTcells.means[,1]

# Remove this row from the dataset.
gdTcells.means <- gdTcells.means[,-1]

#Create query for gene id >- gene name conversion
qs <- ""
for (i in 1:length(genes)){
  qs[i] <- paste("EnsemblGeneID==", genes[i], sep="")
}
q=paste(qs, collapse=",")

#Convert genes
conversion.table <- molgenis.get("lncrna_GeneInfo", q=q)
#row.names(gdTcells.means) = as.character(conversion.table[do.call(c,lapply(row.names(gdTcells.means), function(x) which(as.character(conversion.table[,1]) == x))),2])

#Set colors and breaks
hmcol = colorRampPalette(brewer.pal(9, "GnBu"))
b <- c(-1.5,0,1.5)
colors <- c("#FFFFFF", "#EFCC8C", "#D7A069", "#C1755B", "#A63A4A")
#png("${outputFile}", width=700, height=700)

#Create heatmap
heatmap.2(as.matrix(gdTcells.means),col=colors, trace="none",margin=c(10,10),scale="none",dendrogram="none", Rowv="NULL", 
          Colv="NULL",key=F, cexRow=1.2, labRow = conversion.table[,2],
          cexCol=1.2, lhei=c(1,8),breaks=c(-2,-0.1,1,5,10,max(gdTcells.means)), colsep=c(1:28), sepcolor="snow2", sepwidth=c(0.01,0.01)
          )
image.plot(matrix(-3:12,16,8),legend.only=TRUE, horizontal=TRUE, col=colors, legend.shrink=0.25, legend.width=0.5,
           axis.args=list(at=c(-3,0,3,6,9,12),labels=c("NA",0,3,5,10,"MAX")),smallplot=c(0.51,0.9,0.96,0.99))
