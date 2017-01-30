# Library: MatrixEQTL 2.1.1

## CRAN
# Check if all libraries are installed, if not install them.
if ( length ( setdiff ( "MatrixEQTL", rownames( installed.packages() ) ) ) > 0 ) {
  install.packages( setdiff ( "MatrixEQTL", rownames ( installed.packages() ) ),
                    repos = "http://cran.xl-mirror.nl")
} else {
  ## Libraries
  # Load Matrix eQTL library.
  library("MatrixEQTL")
  
}

## ctQTL.
# Maps cis- and trans-eQTLs
ctQTL <- function () {
  
  # Set gene and snp position files.
  snpspos = read.table("~/Dropbox/Erik Schutte Internship 2016/Data/eQTL-data/snp.txt", header = TRUE, stringsAsFactors = FALSE)
  # genepos = read.table("~/Dropbox/Erik Schutte Internship 2016/Data/eQTL-data/gencode_genes.txt", stringsAsFactors = FALSE)
  genepos = read.table("~/Dropbox/Erik Schutte Internship 2016/Data/eQTL-data/genepos.txt", stringsAsFactors = FALSE, header=T)
  snpspos <- snpspos[-which(snpspos$SNP %in% snps.discarded),]
  # Only associations significant at this level will be saved.
  pvOutputThreshold_tra = 1e-5
  pvOutputThreshold_cis = 0.05
  
  # Distance for local gene-SNP pairs
  cisDist = 5e5 # 500 change this

  ### Prepare matrix eqtl
  ## Settings
  # Set the model used for eqtl mapping.
  useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS
  useCase = "timepoints" # timepoints, all
  
  # Change covariates to character() for no covariates.
  if (useCase == 'all') {
    covariates_file_name = "~/Dropbox/Erik Schutte Internship 2016/Data/eQTL-data/covariates_leiden_samples.txt"
  } else {
    covariates_file_name = character()
  }
  
  
  # Error covariance matrix.
  # Set to numeric() for identity.
  errorCovariance = numeric();
  
  # Empty variable for cis- trans- qtls.
  transqtls <- NULL
  cisqtls <- NULL
  print(dim(snps.t))
  if (useCase == 'timepoints') {
    # Load genotype data.
    snps.sd = SlicedData$new()
    snps.sd$CreateFromMatrix(snps.t)
    snps.sd$fileDelimiter = "\t";      # the TAB character
    snps.sd$fileOmitCharacters = "NA"; # denote missing values;
    snps.sd$fileSkipRows = 1;          # one row of column labels
    snps.sd$fileSkipColumns = 1;       # one column of row labels
    snps.sd$fileSliceSize = 2000;      # read file in slices of 2,000 rows
    
    # Analysis genotype data versus gene expression data requries a loop through the 
    # different timepoints measured in the gene expression data.
    time_intervals <- list(t0 = seq(1,dim(GE)[2],4),
                           t10 = seq(2,dim(GE)[2],4),
                           t180 = seq(3,dim(GE)[2],4),
                           t30 = seq(4,dim(GE)[2],4))
    
    # Change gene expression data to 4 segments for the time intervals.
    for (interval in time_intervals) {
      
      # Pattern for output file name.
      pattern_name.cis <- paste(mainDir,subDir,"lncrna_eqtl_cis_",names(time_intervals[interval])[1],".csv",sep="")
      pattern_name.trans <- paste(mainDir,subDir,"lncrna_eqtl_trans_",names(time_intervals[interval])[1],".csv",sep="")
      
      mapping(pattern_name.cis, pattern_name.trans, useCase, snps.sd,
              snpspos, genepos, pvOutputThreshold_cis, pvOutputThreshold_tra, cisDist,
              useModel, covariates_file_name, errorCovariance, time_intervals, interval)
    }
  } else {
    
    # Multiply genotypes from 22 to 88 samples.
    snps.t <- cbind(snps.t, snps.t) # 44
    snps.t <- cbind(snps.t, snps.t) # 88
    
    # Load genotype data.
    snps.sd = SlicedData$new()
    snps.sd$CreateFromMatrix(snps.t)
    snps.sd$fileDelimiter = "\t";      # the TAB character
    snps.sd$fileOmitCharacters = "NA"; # denote missing values;
    snps.sd$fileSkipRows = 1;          # one row of column labels
    snps.sd$fileSkipColumns = 1;       # one column of row labels
    snps.sd$fileSliceSize = 2000;      # read file in slices of 2,000 rows
    
    # Set pattern names.
    pattern_name.cis <- paste(mainDir,subDir,"lncrna_interactive_eqtl_cis_all_timepoints.csv",sep="")
    pattern_name.trans <- paste(mainDir,subDir,"lncrna_interactive_eqtl_trans_all_timepoints.csv",sep="")
    
    # Do the mapping.
    mapping(pattern_name.cis, pattern_name.trans, useCase, snps.sd,
            snpspos, genepos, pvOutputThreshold_cis, pvOutputThreshold_tra,
            cisDist, useModel, covariates_file_name, errorCovariance)
  }
  
  
  
  
  #transqtls <- format(transqtls,scientific=TRUE)
  #cisqtls <- format(cisqtls,scientific=TRUE)
  
  cat("Writing to file...\n\n")
  #write.csv( transqtls,file=paste(file.path(mainDir, subDir, "Trans", fsep = ""), "/lncrna_eqtls_trans.csv", sep="") )
  #write.csv( cisqtls,file=paste(file.path(mainDir, subDir, "Cis", fsep = ""), "/lncrna_eqtls_cis.csv", sep="") )
  write.csv( MAFGenotype, file=paste(file.path(mainDir, subDir, fsep = ""), "/MAFGenotype.csv", sep="")  )
  
  cat("Wrote file(s) to: \n",paste(file.path(mainDir, subDir, "Trans", fsep=""), sep=""),
      paste(file.path(mainDir, subDir, "Cis", fsep="") ) )
  
  cat("\n\nFINISHED\n\n")
}

### Function
## eQTL mapping.
# Maps using MatrixEQTL.
mapping <- function(pattern_name.cis=pattern_name.cis, pattern_name.trans=pattern_name.trans,
                    useCase, snps.sd, snpspos, genepos, pvOutputThreshold_cis, pvOutputThreshold_tra,
                    cisDist, useModel, covariates_file_name, errorCovariance, time_intervals = NULL, interval = NULL ) {
  # Create output files.
  file.create(pattern_name.cis)
  file.create(pattern_name.trans)
  
  # Output file name and location.
  output_file_name_cis = pattern_name.cis
  output_file_name_tra = pattern_name.trans
  
  # Load gene expression data .
  gene = SlicedData$new();
  if ( useCase == 'timepoints') {
    gene$CreateFromMatrix(GE[,interval])
  } else {
    gene$CreateFromMatrix(GE)
  }
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
  if( useCase == 'all') {
    if(length(covariates_file_name)>0) {
      cvrt$LoadFile(covariates_file_name);
    }
  }
  
  ## Run the analysis.
  me = Matrix_eQTL_main(
    snps = snps.sd,
    gene = gene, 
    cvrt = cvrt,
    output_file_name     = output_file_name_tra,
    pvOutputThreshold     = pvOutputThreshold_tra,
    useModel = useModel, 
    errorCovariance = errorCovariance, 
    verbose = TRUE, 
    output_file_name.cis = output_file_name_cis,
    pvOutputThreshold.cis = pvOutputThreshold_cis,
    snpspos = snpspos, 
    genepos = genepos,
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
  trans <- me$trans$eqtls
  cis <- me$cis$eqtls
  
  if (useCase == 'timepoints') {
    # Also save as Rdata file for Trans eQTLs.
    destination=paste("~/Dropbox/Erik Schutte Internship 2016/Results/Trans/GENCODE/",names(time_intervals[interval])[1],"/gsTcell_",names(time_intervals[interval])[1],"_trans-eQTLs_0.05.RData",sep="")
    save(trans, file=destination)
    
    # Also save as Rdata file for cis eQTLs.
    destination=paste("~/Dropbox/Erik Schutte Internship 2016/Results/Cis/GENCODE/",names(time_intervals[interval])[1],"/gsTcell_",names(time_intervals[interval])[1],"_cis-eQTLs_0.05.RData",sep="")
    save(cis, file=destination)
    
  } else {
    destination=paste("~/Dropbox/Erik Schutte Internship 2016/Results/Cis/GENCODE/gsTcell_all_timepoints_interactive-cis-eQTLs_0.05.RData",sep="")
    save(cis, file=destination)
    
    destination=paste("~/Dropbox/Erik Schutte Internship 2016/Results/Cis/GENCODE/gsTcell_all_timepoints_interactive-trans-eQTLs_0.05.RData",sep="")
    save(trans, file=destination)
  }
}