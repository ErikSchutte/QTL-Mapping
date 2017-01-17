##########################################################
#Authors: Yang Li & Marije van der Geest & Erik Schutte  #
#                                                        #
#Generates heatmap for gluten specific T-cell expression #
#per gene under different conditions (timpoints 0, 10,   #
#30, 180)                                                #
##########################################################

########
# Libs #
########
library("gplots")
library("RColorBrewer")


#Create connection with molgenis
source("http://localhost:8080/molgenis.R")
molgenis.login("admin", "admin")

#Get and process genes (input)
id.my <- "ENSG00000131748,ENSG00000173991,ENSG00000141744,ENSG00000161395,ENSG00000141736,ENSG00000265178,ENSG00000264198,ENSG00000264663,ENSG00000226117,ENSG00000264968,ENSG00000265799,ENSG00000238793,ENSG00000221555,ENSG00000141741,ENSG00000141738,ENSG00000161405,ENSG00000186075,ENSG00000073605,ENSG00000172057,ENSG00000204913,ENSG00000167914,ENSG00000108344,ENSG00000108342,ENSG00000008838,ENSG00000126351,ENSG00000126368,ENSG00000188895,ENSG00000108349,"
genes <- unlist(strsplit(id.my, '[,]'))

#Create query for filtering data
qs <- ""
for (i in 1:length(genes)){
  qs[i] <- paste("EnsemblGeneID==", genes[i], sep="")
}
q=paste(qs, collapse=",")

#Get data
groupMean <- molgenis.get("lncrna_data_Tcell123samples", q=q)

#Format data
rownames(groupMean) <- groupMean[,1]
groupMean <- groupMean[,-1]
data.plot <- groupMean[genes,]

#Remove duplicates
genes.remove <- which(row.names(data.plot) != genes)
if (length(genes.remove) != 0){
  genes <- genes[-genes.remove]
  data.plot <- na.omit(data.plot)
}

#Create query for gene id -> gene name conversion
qs <- ""
for (i in 1:length(genes)){
  qs[i] <- paste("EnsemblGeneID==", genes[i], sep="")
}
q=paste(qs, collapse=",")
#Convert genes
conversion.table <- molgenis.get("lncrna_GeneInfo", q=q)
#row.names(data.plot) = as.character(conversion.table[do.call(c,lapply(row.names(data.plot), function(x) which(as.character(conversion.table[,1]) == x))),2])
#row.names(data.plot) = conversion.table[,1]

#Set colors
hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(100)
#Create image
#png("Rplot.png",width=500, height=400)
seq_cols = seq(0,length(colnames(data.plot)),length(colnames(data.plot))/4)
seq_names = seq( 0, length(colnames(data.plot)), (length(colnames(data.plot))/4.5) )
col.names <- c(rep("",length(colnames(data.plot))))

col.names[seq_names]<-c("t=0","t=10","t=30","t=180")
heatmap.2(as.matrix(data.plot), col = hmcol, trace="none", rowsep=1:nrow(data.plot),
          scale="row",dendrogram="row",Colv=FALSE,colsep=seq_cols,srtCol=0,
          sepcolor="white", cexRow=1, cexCol=1,labCol=col.names,
          labRow = conversion.table[,2], sepwidth=c(0.1,0.00001) )
