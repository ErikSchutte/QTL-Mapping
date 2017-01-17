##########################################################
#Authors: Marije van der Geest & Erik Schutte            #
#                                                        #
#Generates heatmap for rpkm expression data for          #
#different cell types                                    #
##########################################################

########
# Libs #
########
library("gplots")
library("fields")

## Set-up
# Make connection with Molgenis.
source("http://localhost:8080/molgenis.R")

# Verbose - Sample login username "admin", password "admin" for localhost.
molgenis.login("admin", "admin")

## Get gene information
# Prepare query for name conversion of input genes.
# Verbose - Sample input "ENSG00000115598,ENSG00000115602,ENSG00000115604,ENSG00000205716,"
genes <- strsplit("ENSG00000115598,ENSG00000115602,ENSG00000115604,ENSG00000205716,", '[,]')
query <- lapply(genes, FUN=function(geneName) paste("EnsemblGeneID==",geneName, collapse=" or ", sep=""))

# Convert names from 'EnsemblGeneID' to 'AssociatedGeneName'.
name.conversion <- molgenis.get("lncrna_GeneInfo", q=query)

# Select the Gene names.
genes <- as.character(name.conversion[,2])

## Get expression data
# Prepare query for expression data.
query <- lapply(genes, FUN=function(geneName) paste("AssociatedGeneName==",geneName, collapse=" or ", sep=""))

# Get expression data from selected genes.
dane <- molgenis.get("lncrna_data_rpkm7CT", q=query)

## Filter data
# Check which rownames are actual expression data, and which are 'column names',
# this is because the query adds all the results together as one big matrix.
prep.rownames.index <- which(dane[,1]!="AssociatedGeneName")

# Only select expression data.
dane <- dane[prep.rownames.index,]

# Store the 'AssociatedGeneName' for row names.
rownames.dane <- dane[,1]

# Remove the 'AssociatedGeneName' column.
dane <- dane[,2:8]

# Create a separator for the rows.
mysep <- c(1:nrow(dane))              

# Remove factor levels.
dane <- as.matrix(dane)
class(dane) <- "double"

## Create image
# Set colors.
colors <- c("#FFFFFF", "#EFCC8C", "#D7A069", "#C1755B", "#A63A4A")

# Create temporary output file.
#png("${outputFile}", width=500, height=500)

# Plot the heatmap.
heatmap.2(as.matrix(dane),col=colors, trace="none",margin=c(10,10),scale="none",dendrogram="none", Rowv="NULL", 
          Colv="NULL",key=F, cexRow=1.2,
          cexCol=1.2, lhei=c(1,8), rowsep=mysep,breaks=c(-2,-0.1,1,5,10,max(dane)), colsep=c(1:28), sepcolor="snow2", sepwidth=c(0.01,0.01),
          labRow = rownames.dane, labCol=c("NK cells","B-cells","Monocytes","T-memory cells","CD4 T-cells","CD8 T-cells","Granulocytes"))

# Plot legend.
image.plot(matrix(-3:12,16,8),legend.only=TRUE, horizontal=TRUE, col=colors, legend.shrink=0.25, legend.width=0.5,
           axis.args=list(at=c(-3,0,3,6,9,12),labels=c("NA",0,3,5,10,"MAX")),smallplot=c(0.51,0.9,0.96,0.99))