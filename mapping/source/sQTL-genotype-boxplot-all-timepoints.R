library(reshape2)
library(ggplot2)
# Load sqtls mapped with Matrix eQTL.
sqtls <- read.table("~/Dropbox/Erik Schutte Internship 2016/Data/Leafcutter/intron_splice_sqtl_mapping.csv",
           header=T, stringsAsFactors = F)

load("~/Dropbox/Erik Schutte Internship 2016/Data/Leafcutter/ordered_sqtl_data.Rdata")

source("~/Dropbox/Erik Schutte Internship 2016/Code/QTL-Mapping/mapping/source/prepare_df_sqtls.R")
sqtls.ma <- lapply(list(sqtls), prepare_df_sqtls)[[1]]

sqtls.tmp <- sqtls.ma[1:15,]
ll <- apply(sqtls.tmp, 1, function(sqtl) {
    # Set each row as a dataframe instead of a vector
    sqtl <- data.frame(t(sqtl))
    # sqtl <- sqtls.ma[1,]
    
    # Position current gene in the expression file.
    cluster_count <- introns.filtered[ which( sqtl$gene == rownames( introns.filtered ) ), ]
    
    # Add a column to the gene expression dataframe and fill it with timepoints for later use.
    df.cc <- data.frame(t(cluster_count))
    df.cc <- cbind.data.frame(df.cc, time=c(0))
    df.cc[times$t0,2] <- "t0"
    df.cc[times$t10,2] <- "t10"
    df.cc[times$t30,2] <- "t30"
    df.cc[times$t180,2] <- "t180"
    df.cc.melt <- melt(df.cc)
    
    # Set sample names.
    df.cc.melt[,2] <- colnames(cluster_count)
    
    # Run an index on which snp corrosponds to which genotypes column.
    index = sapply(sqtl$snps, function(snp) {
      which(snp == colnames(snps))
    })
 
    # Retrieve all genotypes from all samples for the specified snps.
    genotypes <- snps.geno[-1,index] # because of missing sample now 1:21

    # Transpose and to data frame.
    genotypes <- t(data.frame(genotypes))

    # Multiply the genotypes 4 times.
    genotypes <- cbind.data.frame(genotypes, genotypes, row.names = NULL) # 44x / 42x
    genotypes <- cbind.data.frame(genotypes, genotypes, row.names = NULL) # 88x / 84x

    # Melt the transposed genotypes.
    genotypes.melt <- melt(t(genotypes))
    genotypes.melt <- genotypes.melt[ order(as.character(genotypes.melt[,1])), ]

    # Bind the dataframes together.
    df.melt <- cbind.data.frame(timepoints=df.cc.melt[,1],
                                expression=df.cc.melt[,3],
                                sample=df.cc.melt[,2],
                                genotypes=genotypes.melt)
  
    # Create an order for the timepoints.
    timepoints.order <- c("t0","t10","t30","t180")
    
    # Order the factor levels for the timepoints.
    df.melt$timepoints <- factor(df.melt$timepoints, levels=timepoints.order)
    
    # Create an order for the genotypes.
    genotype.levels <- levels(df.melt$genotypes.value)
    major.allele = strsplit(as.character(sqtl$minor_major),"_")[[1]][2]
    genotype.order <- c()
    if (length(genotype.levels) == 2) { # Only 2 levels for genotypes.
      
      one <- strsplit(genotype.levels,"")[[1]] # first level
      two <- strsplit(genotype.levels,"")[[2]] # second level
      
      if (major.allele == one[1] && major.allele == one[2]) { # first level is homozygous major allele.
        print("one is major")
        major = genotype.levels[1]
      } else if (major.allele == one[1] && major.allele != one[2] | major.allele != one[1] && major.allele == one[2]) { # first level is heteryzygous
        print("one is hetero")
        heter = genotype.levels[1]
        genotype.order <- c(genotype.order, one)
      } else {
        print("WARNING: Third genotype should not exists here.") 
      }
      
      if (major.allele == two[1] && major.allele == two[2]) { # first level is homozygous major allele.
        print("two is major")
        major = genotype.levels[2]
      } else if (major.allele == two[1] && major.allele != two[2] | major.allele != two[1] && major.allele == two[2]) { # first level is heteryzygous
        print("two is heter")
        heter = genotype.levels[2]
      } else {
        print("WARNING: Third genotype should not exists here.") 
      }
      
      genotype.order <- c(major, heter)
      
    } else { # 3 levels for genotypes.
      one <- strsplit(genotype.levels,"")[[1]] # first level
      two <- strsplit(genotype.levels,"")[[2]] # second level
      three <- strsplit(genotype.levels,"")[[3]] #third level
      if (major.allele == one[1] && major.allele == one[2]) { # first level is major allele.
        print("one is major")
        major = genotype.levels[1]
      } else if (major.allele == one[1] && major.allele != one[2] | major.allele != one[1] && major.allele == one[2]) { # first level is heteryzygous
        print("one is hetero")
        heter = genotype.levels[1]
        genotype.order <- c(genotype.order, one)
      } else {
        print("one is minor")
        minor = genotype.levels[1]
      }
      
      if (major.allele == two[1] && major.allele == two[2]) { # first level is homozygous major allele.
        print("two is major")
        major = genotype.levels[2]
      } else if (major.allele == two[1] && major.allele != two[2] | major.allele != two[1] && major.allele == two[2]) { # first level is heteryzygous
        print("two is heter")
        heter = genotype.levels[2]
      } else {
        minor = genotype.levels[2]
      }
      
      if (major.allele == three[1] && major.allele == three[2]) { # first level is homozygous major allele.
        print("three is major")
        major = genotype.levels[3]
      } else if (major.allele == three[1] && major.allele != three[2] | major.allele != three[1] && major.allele == three[2]) { # first level is heteryzygous
        print("three is heter")
        heter = genotype.levels[3]
      } else {
        minor = genotype.levels[3]
      }
      
      genotype.order <- c(major, heter, minor)
    }

    # Order the factor levels for the genotypes.
    df.melt$genotypes.value <- factor(df.melt$genotypes.value, levels=genotype.order)
    
    # Order data on factor levels.
    df.melt <- df.melt[order(df.melt$timepoints),]
    df.melt <- df.melt[order(df.melt$sample),]
    
    # Save ggplot in variable and plot.
    mi <- min( df.melt$expression ) - 1
    ma <- max(df.melt$expression ) + 1
    p <- ggplot(data=df.melt, aes(x=genotypes.value, y=expression, group=genotypes.value) ) +
      geom_boxplot(aes( fill=genotypes.value), outlier.shape=NA ) +
      geom_point( position=position_jitter(width=0.15),colour = "darkgrey") +
      coord_cartesian( ylim = c( mi,ma ) ) +
      ggtitle( paste(sqtl$gene, " - ", sqtl$snps, sep = ""),
               subtitle = paste( "P-value: ", sqtl$pvalue ) ) +
      theme( plot.title = element_text( size = rel(1.6), hjust = 0.5 ),
             plot.subtitle = element_text(size = rel(1), hjust = 0.5 ) ) +
      xlab(paste("Genotypes",sep="")) + ylab("Scaled+Norm. Intron Count") +
      scale_fill_discrete( name="Genotypes",
                           labels=paste( names( table( df.melt$genotypes.value ) ),"(", table( df.melt$genotypes.value )/4, ")", sep ="") )
    p  + facet_wrap( ~ timepoints, scales="free")
    
    cluster=gsub(as.character(sqtl$gene), pattern=":", replacement="_")
    ggsave(filename=paste("~/Dropbox/Erik\ Schutte\ Internship\ 2016/Images/sQTLs/",cluster,"_",sqtl$snps,"_",sqtl$minor_major,".pdf",sep=""),
           plot=last_plot(), device = "pdf")

    
    pdf(paste("~/Dropbox/Erik\ Schutte\ Internship\ 2016/Images/sQTLs/",cluster,"_",sqtl$snps,"_",sqtl$minor_major,"_","correlation.pdf",sep=""),
        width=7, height=5)
    par(mfrow=c(2,2))
    geno <- data.frame(snps[,which(colnames(snps) == sqtl$snps)])
    timepoints <- c("t0", "t10", "t30", "t180")
    for ( time in timepoints ) {
      names1 <- as.character(df.melt[which(df.melt$timepoints == time),3])
      plot( geno[names1,], df.melt[which(df.melt$timepoints == time),2], ylim=c(mi,ma) , main = paste("Gene: ", sqtl$gene, "\nSNP: ", sqtl$snps ),
            xlab = "Genotype", ylab = "Gene expression", sub = time, type = "p", pch=20)
      abline(lm(df.melt[which(df.melt$timepoints == time),2] ~ geno[names1,] ), col="red" )
      abline(h = 0, col="grey")
    }
    dev.off()
    
    # Set covariates.
    age = unlist(covariates[1,])
    gender = unlist(covariates[2,])
    time.interval = unlist(covariates[3,])
    
    # Calculate residuals.
    resi=lm(unlist(df.melt$expression)~age+gender+time.interval)$residuals
    names(resi) <- df.melt$sample
    
    # Get genotypes.
    gt1=geno[names(resi),]
    
    # Create matrix for residuals and genotypes.
    data1=cbind(resi, gt1)
    print(cor.test(data1[,2], data1[,1])$p.value)

})

