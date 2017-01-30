## Library
library(gtools)

## Data preparation.
data_prep <- function () {
  ### Data-preparation
  ## Transform the snps.
  # Verbose.
  if ( VERBOSE == TRUE ) {
    # Verbose - Print transforming SNPs.
    cat("Transforming SNPs\n\n")
  }
  
  ##### Load input data ####
  load("~/Dropbox/Erik Schutte Internship 2016/Data/eQTL-data/clones_CeD_genotypes_noDuplicates.Rdata")
  load("~/Dropbox/Erik Schutte Internship 2016/Data/eQTL-data/expr_vst_condition_patient_rmBatch_88samples.Rdata")
  
  ##### Transform names ####
  ## Set data to shorter values.
  ge <- expr.vst.rmBatch
  snps <- mergedGenotypes_clones2
  
  ## For readability change the matrix name and remove Norway data.
  # Norway data is excluded from our eqtl analyis.
  snps <- snps[-(23:27),]
  
  ## Add a sample that hasn't been added yet.
  # There is a sample that is not included in the genotype data but has the same genotype as TCC-01-2 so we add that.
  tmp.df <- t(data.frame("TCC-01-3"=snps[15,]))
  snps <- rbind(snps, tmp.df)
  
  # Set the rowname for the new sample to TCC-01-3 instead of TCC.01.3.
  rownames(snps)[which(rownames(snps) == "TCC.01.3")] <- "TCC-01-3"
  
  ## Order Names.
  # Order the sample names.
  snps <- snps[order(rownames(snps)),]
  
  # Add a row to the data frame with MAF where the alleles will be stored.
  tmp.df <- t(data.frame(MAF=rep("",length(colnames(snps)))))
  snps.genotype <- snps # Copy of snps that won't be changed to frequencies.
  
  # Change the sample names to the corrosponding format for the Gene Expression file.
  if ( length(grep("\\.",colnames((ge))) != 0) ) {
    colnames(ge) <- gsub("\\.", "_", colnames(ge))
  } else if ( length(grep("\\-",colnames((ge))) != 0) ) {
    colnames(ge) <- gsub("\\-", "_", colnames(ge))
  } 
  # Change column names to match sample names.
  colnames(ge) <- gsub("batch[0-9]+__", "", colnames(ge))
  ge <- ge[,mixedsort(colnames(ge))]
  colnames(ge) <- gsub("__t[0-9]+", "", colnames(ge))
  
  # Change the sample names to the corrosponding format for the Genotype file.
  if ( length(grep("-",rownames(snps)) != 0) ) {
    rownames(snps) <- gsub("\\-", "_", rownames(snps))
    rownames(snps.genotype) <- gsub("\\-", "_", rownames(snps.genotype))
  }
  
  # Remove those samples that don't occur in the gene expression data. 
  # “TCC-09-1", ”TCC-08-1”, "TCC-21-01" and "TCC-20-01".
  com.samples <- intersect(rownames(snps), colnames(ge))
  
  ## Re-arange data frame.
  # If we re-arrange without saving the ordering for the samples we lose it by using the indices.
  sample.order <- rownames(snps) # sample order
  snps <- snps[com.samples,] 
  snps.genotype <- snps.genotype[com.samples,]
  snps.genotype <- rbind(snps.genotype, tmp.df) # after re-arangement bind alleles to genotypes.
  
  ## Order the gene expression data
  ge <- ge[,which(colnames(ge)%in%com.samples)]
  
  ##### Preparation ####
  ## Convert genotypes to 0, 1 and 2. 0 ref/ref, 1 ref/mut and 2 mut/mut.
  # Keeps track of the snps that are causing noise in the data and are not representable.
  snps.discarded <- c()
  snps.discarded.pos <- c()
  snps.discarded.counter <- 0
  # For every column in the genotype data.
  for (i in 1:ncol(snps)) {
    
    # Verbose
    if ( VERBOSE == TRUE ) {
      # Verbose - Print current column and column name.
      cat("#",i,"- current column: ",colnames(snps)[i],"\n")
    }
    # Save genotype and genotype count for each column.
    alleles.all <- table(snps[,i])
    
    # A minimum of 3 samples for a genotype is required to effectively calculate correlation.
    # Therefore we are removing the rows where this is not the case.
    if ( length(alleles.all) < 3 & min(alleles.all) < 2) {
      
      # Set an index for the SNP, so that it can be removed later.
      snps.discarded.pos <- c(snps.discarded.pos,i)
      
      # Store the discarded SNP for later review.
      snps.discarded <- c(snps.discarded, colnames(snps)[i])
      snps.discarded.counter <- snps.discarded.counter + 1
    }
    
    # List with all the major, minor, and heterozygote genotypes and their respecitve counts.
    allele.list <- list(major <- list(genotype =NULL,count=0),
                        minor <- list(genotype =NULL,count=0),
                        hetero <- list(genotype =c(),count=c()),
                        other <- list(genotype =NULL,count=0))
    
    # Temp value for checking which genotype occurs the most.
    highest_count_homozygote <- 0
    
    # For every allele that the genotype is counted for.
    for (a in 1:length(alleles.all)) {
      
      # Extract the genotype.
      geno <- names(alleles.all[a])
      
      # Verbose
      if ( VERBOSE == TRUE ) {
        # Verbose - Print current allele.
        cat("Current allele = ", geno,'\n')
      }
      
      # Extract the number of occurrences.
      count <- alleles.all[[a]]
      
      # Split the genotype 
      alleles <- strsplit(geno, "")
      
      #Check if the first allele is the same as the second
      # thus checking for heterozygote or homozygote genotype.
      if (identical(alleles[[1]][1],alleles[[1]][2]) & geno != '00') {
        
        # Verbose
        if ( VERBOSE == TRUE ) {
          # Verbose - Print if alleles are identical, and thus are homozygous.
          cat('Alleles are identical, sorted as homozygote\n')
        }
        
        # If the highest_count_homozygote == 0 (the first iteration).
        if (highest_count_homozygote == 0) {
          
          # Verbose
          if ( VERBOSE == TRUE ) {
            # Verbose - Print which genotype and count is added as major allele.
            cat("Adding ",geno," with count ", count, "as major allele\n")
          }
          
          # For the first iteration, the first genotype's occurences are the most occured
          # genotype and will be added on the first position.
          highest_count_homozygote <- count
          allele.list[[1]]$genotype <- geno
          allele.list[[1]]$count <- count
          
          # If the current genotype has more occurences than the lastly noted genotype,
          # it is automatically the major allele on the first position.  
        } else if (highest_count_homozygote < count) {
          
          # Verbose
          if ( VERBOSE == TRUE ) {
            # Verbose - Print which genotype and count is added as new major allele,
            # and which genotype and count previously occupied that spot.
            cat("Adding ",geno," with count ", count, "as major allele\nPreviously ",
                allele.list[[1]]$genotype," with count ", allele.list[[1]]$count,'\n')
          }
          
          # Altering the position from the previous current major allele to minor allele
          # in the allele.list.
          allele.list[[2]]$genotype <- allele.list[[1]]$genotype
          allele.list[[2]]$count <- allele.list[[1]]$count
          
          # Setting the new major allele on the first position.
          highest_count_homozygote <- count
          allele.list[[1]]$genotype <- geno
          allele.list[[1]]$count <- count
          
        } else {
          # If the highest occuring genotype already found has more occurences,
          # the resulting genotype will be the minor allele which is the 2nd position
          # in the list.
          allele.list[[2]]$genotype <- geno
          allele.list[[2]]$count <- count
          
          # Verbose
          if ( VERBOSE == TRUE ) {
            # Verbose - Print which genotype and count is added as minor allele.
            cat("Adding ",geno," with count ", count, "as minor allele\n")
          }
        }
        
      } else {
        if (identical(geno,'00')) {
          # Separating the '00'/NA's from the data, these will be represented as 9.
          allele.list[[4]]$genotype <- geno
          allele.list[[4]]$count <- count
          
        } else {
          # If the genotype is not '00' it is automatically a heterozygote.
          allele.list[[3]]$genotype <- c(allele.list[[3]]$genotype, geno)
          allele.list[[3]]$count <- c(allele.list[[3]]$count, count)
          
          # Verbose
          if ( VERBOSE == TRUE ) {
            # Verbose - Print if genes are not identical, and thus are heterozygote.
            cat('alleles are not identical, heterozygote\n')
            # Verbose - Print which genotype and count is added as heterozygous.
            cat("Adding ",geno," with count ", count, "as heterozygote allele\n")
          }
        }
      }
    }
    
    # For every row in the genotype data.
    for (j in 1:nrow(snps.genotype) ) {
      
      # First we want to create a dict that takes the inverse of the MAJOR alleles we put there,
      # since we do not save all the minor alleles.
      #list()
      
      # If J is smaller than the maximum lengt of rows, e.g. all the samples minus the last row
      # where the MAF will be stored.
      if (j < (nrow(snps.genotype) ) ) {
        cat("Current is: ",j,"\n")
        # Check the current snp agains the determined major, minor and heterozygote genotypes,
        # and determine the allele frequency.
        snp.freq <- change_allele_to_frequencies(snps[j,i],allele.list)
        
        # Verbose.
        if ( VERBOSE == TRUE ) {
          # Verbose - Print the current snip' genotype and the corrosponding frequency.
          cat("Changing current snip: ",snps[j,i]," to: ", snp.freq,"\n")
        }
        if ( length(allele.list[[3]]$genotype) > 1) {
          
          if ( snps[j,i] != allele.list[[3]]$genotype[1] && snps[j,i] == allele.list[[3]]$genotype[2]) {
            snps.genotype[j,i] <- allele.list[[3]]$genotype[1]
          }
        }
        
        
        # Save the SNP frequency at previous occupied genotype location.
        snps[j,i] <- snp.freq
      }
      # J is equal to the length of rows for the data frame snps. meaning it is hte last row that will contain
      # the MAF.
      else {
        cat("Other is: ",j,"\n")
        # We don't always have the minor allele, so instead we look at the major and heterozygous genotypes.
        # E.g. all 22 samples have the following genotyeps for snp X, TT (19x) CT (3x) CC(0x).
        # Because we can't count CC we have no genotype CC and thus have to 'create' it.
        major.allele = strsplit(allele.list[[1]]$genotype,"")[[1]] # Contains both alleles, e.g. T and T.
        
        if ( length(allele.list[[3]]$genotype) < 1 && length(allele.list[[2]]$genotype) < 1 ) {
          
          # In the off case that there are TT(22x) CT(0x) CC(0x)
          minor.allele <- "?"
        } else if ( length(allele.list[[3]]$genotype) < 1 && length(allele.list[[2]]$genotype) > 0 ) {
          
          # In the off case that there are TT(17x) CT(0x) CC(5x)
          minor <- allele.list[[2]]$genotype
          minor.alleles <- strsplit(minor,"")[[1]]
          minor.allele <- minor.alleles[1]
          
        } else {
          
          heter.allele = strsplit(allele.list[[3]]$genotype,"")[[1]] # Contains both alleles, e.g. C and T.
          
          if ( heter.allele[1] == major.allele[1] ) {
            
            # If the first allele, in our example C matches the major allele T we know that he minor allele must be
            # the second allele of the heterozygous alleles.
            minor.allele <- heter.allele[2]
            
          } else if ( heter.allele[2] == major.allele[1] ) {
            
            # If the second allele, in our exmaple T matches the major allele T we know that the minior allele must be 
            # the first allele of the heterozgyous alleles.
            minor.allele <- heter.allele[1]
          }
        }
        
        # Save the genotypes in the data frame as MINOR / MAJOR.
        alleles <- paste(minor.allele, "_", major.allele[1], sep="")
        print(alleles)
        snps.genotype[j,i] <- alleles
      }
    }
  }
  
  # Remove discarded snps from data set.
  snps <- snps[,-snps.discarded.pos]
  snps.genotype <- snps.genotype[,-snps.discarded.pos]
  
  # Convert genotype matrix to numeric values for Matrix eQTL analaysis.
  suppressMessages(
    class(snps) <- "numeric"
  )
  
  # Transfer the genotype matrix so that the columns 'align' with the gene expression matrix.
  snps.t <- t(snps)
  snps.genotype.t <- t(snps.genotype)
}

save(snps.t, file="~/Dropbox/Erik Schutte Internship 2016/Data/eQTL-data/gsTcell_genotype_frequencies_transformed.Rdata")
save(snps.genotype.t, file="~/Dropbox/Erik Schutte Internship 2016/Data/eQTL-data/gsTcell_genotypes_transformed.Rdata")
save(ge, file="~/Dropbox/Erik Schutte Internship 2016/Data/eQTL-data/gsTcell_gene_expr_vst_88samples.Rdata")
