## sQTL.
# Maps splice eQTLs.
sQTL <- function () {
  ### Pre-processing
  ## Load data.
  # Set working directory.
  #setwd("~/Dropbox/Erik\ Schutte\ Internship\ 2016/")
  cat("Starting splice QTL analysis..\n")
  
  # Load in the data files.
  trans.exp.f = "/groups/umcg-wijmenga/tmp04/umcg-eschutte/projects/gluten_specific_Tcells/data/gsTcell_TranscriptExpression.csv"
  cat("Loading gene bed file..\n")
  gene.bed.f = "/groups/umcg-wijmenga/tmp04/umcg-eschutte/projects/gluten_specific_Tcells/data/genes.bed"
  
  ## Order data.
  # Ordered genotype file should be compressed and indexed if not done before.
  cat("Indexing genotype file..\n")
  if ( file.exists("/groups/umcg-wijmenga/tmp04/umcg-eschutte/projects/gluten_specific_Tcells/data/genotype.indexed.f.Rdata") ) {
    load("/groups/umcg-wijmenga/tmp04/umcg-eschutte/projects/gluten_specific_Tcells/data/genotype.indexed.f.Rdata")
  } else{
    
    # Load the genotype Information without the samples, they are generated in this script so we are going to do it on the fly.
    genoTypeInformation.unfinished = read.table("/groups/umcg-wijmenga/tmp04/umcg-eschutte/projects/gluten_specific_Tcells/data/genoTypeInformation.tsv",header=T,sep="\t")
    # There are less snps in the positions file than in my genotype file, find out which are not there and exclude them.
    snps_present = which(genoTypeInformation.unfinished[,4] %in% rownames(snps.t))
    
    # Create a tmp variable to store the new matrix in.
    tmp = genoTypeInformation.unfinished[snps_present,]
    
    # Save an order for the chromosomes.
    chrOrder<-c(paste(1:22,sep=""),"chrX")
    
    # Order the chromosomes factors on the new order.
    tmp[,1] <-factor(tmp[,1], levels=chrOrder)
    
    # Order the data on the new chromosome order.
    tmp <- tmp[order(tmp[,1],tmp[,2]),]
    
    # We can just cbind here, because the rs from snps.t and tmp are orderd the same.
    genoTypeInf <- cbind(tmp, snps.t)
    
    # Remove the Missing ones, or rather replace them with -1.
    genoTypeInf[is.na(genoTypeInf)] <- -1
    
    # Change the factor level so we don't get an NA as chr.
    levels(genoTypeInf[,1])[23] <- 23
    
    # Change the X chromosome to 23.
    genoTypeInf[is.na(genoTypeInf), 1] <- 23
    #colnames(genoTypeInf) <- gsub("\"","",colnames(genoTypeInf))
    
    cat("Loading genotype file..\n")
    # Load the genotype file, we have to do this, because sQTLseeker's logic only accepts a file path instead of a file/matrix or anythings else.
    genotype.f = "/groups/umcg-wijmenga/tmp04/umcg-eschutte/projects/gluten_specific_Tcells/data/genoTypeinf.tsv"
    
    # Write the now finished genotype information matrix to a file, since sqtlseeker only accepts file input.
    write.table( genoTypeInf, file="/groups/umcg-wijmenga/tmp04/umcg-eschutte/projects/gluten_specific_Tcells/data/genoTypeinf.tsv", 
                 sep="\t", row.names = F, quote = F)
    
    # Index the genotype file.
    genotype.indexed.f = index.genotype(genotype.f)
    
    # Save an image of the genotype file.
    save.image(file="/groups/umcg-wijmenga/tmp04/umcg-eschutte/projects/gluten_specific_Tcells/data/genotype.indexed.f.Rdata")
    save(genotype.indexed.f, file="/groups/umcg-wijmenga/tmp04/umcg-eschutte/projects/gluten_specific_Tcells/data/genotype.indexed.f.Rdata")
  }
  
  ## Read data.
  # Import transcript expression, cleaned.
  if ( file.exists("/groups/umcg-wijmenga/tmp04/umcg-eschutte/projects/gluten_specific_Tcells/data/tre.df.Rdata") ) {
    load("/groups/umcg-wijmenga/tmp04/umcg-eschutte/projects/gluten_specific_Tcells/data/tre.df.Rdata")
  } else{
    cat("Reading transcript expression..\n")
    te.df = read.table(trans.exp.f, as.is=T, header=T, sep="\t")
    cat("Prepare transcript expression..\n")
    tre.df = prepare.trans.exp(te.df, min.transcript.exp = 0, min.gene.exp = 0, min.dispersion = 0, verbose = T )
    save.image(file="/groups/umcg-wijmenga/tmp04/umcg-eschutte/projects/gluten_specific_Tcells/data/tre.df.Rdata")
    save(tre.df, file="/groups/umcg-wijmenga/tmp04/umcg-eschutte/projects/gluten_specific_Tcells/data/tre.df.Rdata")
  }
  
  ### Test for gene/SNP associations
  ## Tests.
  # Run test with transcript expression and genotype file and gene coordinates.
  cat("Read gene bed file...\n")
  gene.bed = read.table(gene.bed.f, as.is=T, sep="\t")
  
  # Set column names.
  cat("Set column gene bed file..\n")
  colnames(gene.bed) = c("chr","start","end","geneId")
  
  # Change col names format for tre.df to _ instead of .
  colnames(tre.df) <- sub("\\.","_",colnames(tre.df))
  colnames(tre.df) <- sub("\\.","_",colnames(tre.df))
  
  # Filter Tre.df for samples that don't exist in the genotype file.
  
  ### DYNAMIC SQTL ###
  # Save index for each time point.
  #sort(colnames(tre.df[3:length(colnames(tre.df))]))
  timepoints = list(t0 = grep("_t0",colnames(tre.df)),
                    t10 = grep("_t10",colnames(tre.df)),
                    t30 = grep("_t30",colnames(tre.df)),
                    t180 = grep("_t180",colnames(tre.df)))
  times <- c("t0","t10","t30","t180")
  iteration = 1
  for ( timepoint in timepoints) {
    cat("Do sqtlseeker dynamic..\n")
    
    # For every timepoint create a new data frame, set the first two columns to trId and geneId.
    tre.df.tp <- tre.df[,1:2]
    
    # Add the timepoints according to a grepped index to the new table.
    tre.df.tp <- cbind(tre.df.tp, tre.df[,timepoint])
    
    # Change the column names, so that they match with the sample names from the genotype file.
    colnames(tre.df.tp) <- sub("batch[0-9]+_T","T",colnames(tre.df.tp))
    colnames(tre.df.tp) <- sub("_t[0-9]+","",colnames(tre.df.tp))
    
    # Filter the samples from the transcript expression file that don't exist in the genotype file.
    noGenoTypeSampleIndex <- which(!colnames(tre.df.tp) %in% colnames(snps.t))
    tre.df.tp <- tre.df.tp[,-noGenoTypeSampleIndex[3:length(noGenoTypeSampleIndex)]]
    
    # There is apperantly a sample in batch 2 that is also in batch 4. That only comes to the front if I execute this
    # twice, I don't know why this happens but this solves this.
    noGenoTypeSampleIndex <- which(!colnames(tre.df.tp) %in% colnames(snps.t))
    tre.df.tp <- tre.df.tp[,-noGenoTypeSampleIndex[3:length(noGenoTypeSampleIndex)]]
    
    # Start sqtl.seeker.
    res.df.dynamic = sqtl.seeker.modified(tre.df.tp, genotype.indexed.f, gene.bed, svQTL=F, genic.window=250000, verbose=F)
    
    cat("write to table..\n")
    
    # Write the output to a tsv file.
    write.table(res.df.dynamic,
                file=paste("/groups/umcg-wijmenga/tmp04/umcg-eschutte/projects/gluten_specific_Tcells/R/Splice/sQTLs-dynamic-",times[iteration],"all.tsv",sep=""),
                quote=FALSE, row.names=FALSE, sep="\t")
    
    # Save an image so we can work with the data later if we want to.
    save.image(file=paste("/groups/umcg-wijmenga/tmp04/umcg-eschutte/projects/gluten_specific_Tcells/R/Splice/res.df.dynamic-",times[iteration],"-.Rdata",sep=""))
    save(res.df.dynamic, file=paste("/groups/umcg-wijmenga/tmp04/umcg-eschutte/projects/gluten_specific_Tcells/R/Splice/res.df.dynamic-",times[iteration],"-.Rdata",sep=""))
    
    # Execute sqtls, a function that checks the qvalue of ech sqtl.
    sqtls.df = sqtls(res.df = res.df.dynamic, FDR = 1, out.pdf=paste("/groups/umcg-wijmenga/tmp04/umcg-eschutte/projects/gluten_specific_Tcells/R/Splice/sQTLs-dyanmic-",times[iteration],"-FDR1.pdf", sep="") )
    
    cat("\nProcessed: ",times[iteration])
    
    # Increase iteration time, mainly for name.
    iteration <- iteration + 1
  }
  
  ### STATIC SQTL ###
  tre.df.static <- tre.df
  cat("Do sqtlseeker static..\n")
  colnames(tre.df.static) <- sub("batch[0-9]+_T","T",colnames(tre.df.static))
  colnames(tre.df.static) <- sub("_t[0-9]+","",colnames(tre.df.static))
  res.df.static = sqtl.seeker.modified(tre.df.static, genotype.indexed.f, gene.bed, svQTL=F, verbose=F, genic.window=250000)
  cat("write to table..\n")
  write.table(res.df.static,
              file="/groups/umcg-wijmenga/tmp04/umcg-eschutte/projects/gluten_specific_Tcells/R/Splice/sQTLs-static-all.tsv",
              quote=FALSE, row.names=FALSE, sep="\t")
  sqtls.df = sqtls(res.df.static, FDR = 1, out.pdf="/groups/umcg-wijmenga/tmp04/umcg-eschutte/projects/gluten_specific_Tcells/R/Splice/sQTLs-static-FDR01.pdf")
}