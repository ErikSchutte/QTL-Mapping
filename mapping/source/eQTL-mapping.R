## CRAN
# Check if all libraries are installed, if not install them.
if ( length ( setdiff ( "MatrixEQTL", rownames( installed.packages() ) ) ) > 0 ) {
  supressMessages(install.packages( setdiff ( "MatrixEQTL", rownames ( installed.packages() ) ),
                                    repos = "http://cran.xl-mirror.nl") )
}
## Libraries
# Load Matrix eQTL library.
library("MatrixEQTL")

## eQTL.
# Maps basic eQTLs.
eQTL <- function(){
  
  ### Prepare matrix eqtl
  ## Settings
  # Set the model used for eqtl mapping.
  useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS
  
  # Change covariates to character() for no covariates.
  covariates_file_name = character()
  
  # Only associations significant at this level will be saved.
  pvOutputThreshold = 1e-5;
  
  # Error covariance matrix.
  # Set to numeric() for identity.
  errorCovariance = numeric();
  
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
                         t30 = seq(3,dim(GE)[2],4),
                         t180 = seq(4,dim(GE)[2],4))
  
  # Empty eqtls variable.
  eqtls <- NULL
  
  # Change gene expression data to 4 segments for the time intervals.
  for (interval in 1:length(time_intervals) ) {
    
    # Pattern for output file name.
    pattern_name <- paste(mainDir,subDir,"lncrna_eqtl_basic_",names(time_intervals[interval])[1],".csv",sep="")
    file.create(pattern_name)
    
    # Output file name and location.
    output_file_name = pattern_name
    
    # Load gene expression data.
    gene = SlicedData$new();
    gene$CreateFromMatrix(GE[,time_intervals[[interval]]])
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
    
    ## Run the analysis.
    me = Matrix_eQTL_engine(
      snps = snps.sd,
      gene = gene,
      cvrt = cvrt,
      output_file_name = output_file_name,
      pvOutputThreshold = pvOutputThreshold,
      useModel = useModel, 
      errorCovariance = errorCovariance, 
      verbose = TRUE,
      pvalue.hist = TRUE,
      min.pv.by.genesnp = FALSE,
      noFDRsaveMemory = FALSE);
    
    # Unlink output_file_name (keep it commented to keep eQTL mappings),
    unlink(output_file_name)
    
    # Results.
    cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n\n\n')
    
    # Verbose - Show all detected eQTLs.
    if ( VERBOSE == TRUE ) {
      cat('Detected eQTLs:', '\n')
      show(me$all$eqtls)
    }
    
    eqtls <- rbind(eqtls, me$all$eqtls)
    # Plot the histogram of all p-values.
    #plot(me)
  }
  
  #eqtls.rearanged <- cbind(eqtls[,4],eqtls[,1],eqtls[,2],eqtls[,5])
  #colnames(eqtls.rearanged) <- c("pvalue","SnpRs","EnsembleGeneID","FDR")
  #rownames(eqtls.rearanged) <- rownames(eqtls)
  eqtls <- format(eqtls,scientific=TRUE)
  
  cat("Writing to file...\n\n")
  write.csv(eqtls,file=paste(file.path(mainDir,subDir,"Basic",fsep=""), "/lncrna_eqtls_basic",sep=""))
  
  cat("Wrote file to: ",paste(file.path(mainDir,subDir,"Basic",fsep=""), "/lncrna_eqtls_basic",sep=""),"\n\n")
  cat("FINISHED\n\n")

}
