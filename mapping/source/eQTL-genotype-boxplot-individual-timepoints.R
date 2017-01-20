default.wd <- getwd()
### Data
## Biomart.
# Load biomart data for HNPC Gene names.
biomart <- read.table("~/Dropbox/Erik Schutte Internship 2016/Data/eQTL-data/gencode_genenames.txt",header = F, sep=" ")
colnames(biomart) <- c("gene_id", "gene_name")

## Genotype data
# Load genotype data.
snps <- read.table("~/Dropbox/Erik Schutte Internship 2016/Data/filtered-data/snps.genotype.t.tsv",header=T)
snps <- t(snps)

## Expression data.
# Load expression data.
load("~/Dropbox/Erik Schutte Internship 2016/Data/eQTL-data/gsTcell_vst_88samples.Rdata")
gene.expr <- vst

# Gencode annotated raw counts normalized with VST. However, one sample is missing, remove the entire sample for the time being.
gene.expr <- gene.expr[,-grep("batch\\d+_TCC-17",colnames(gene.expr))]

### PER TIME POINT
## Load data files for eqtl plots.
# Cis eqtls dynamic.
load("~/Dropbox/Erik Schutte Internship 2016/Results/Cis/GENCODE/t0/gsTcell_t0_cis-eQTLs_0.05.RData")
c.0 <- cis[1:5,]
c.0 <- cbind(c.0, data.frame(t.interval=c("_t0"),origin=c("cis")))
load("~/Dropbox/Erik Schutte Internship 2016/Results/Cis/GENCODE/t10/gsTcell_t10_cis-eQTLs_0.05.RData")
c.10 <- cis[1:5,]
c.10 <- cbind(c.10, data.frame(t.interval=c("_t10"),origin=c("cis")))
load("~/Dropbox/Erik Schutte Internship 2016/Results/Cis/GENCODE/t30/gsTcell_t30_cis-eQTLs_0.05.RData")
c.30 <- cis[1:5,]
c.30 <- cbind(c.30, data.frame(t.interval=c("_t30"),origin=c("cis")))
load("~/Dropbox/Erik Schutte Internship 2016/Results/Cis/GENCODE/t180/gsTcell_t180_cis-eQTLs_0.05.RData")
c.180 <- cis[1:5,]
c.180 <- cbind(c.180, data.frame(t.interval=c("_t180"),origin=c("cis")))

# Trans eqtls dynamic.
load("~/Dropbox/Erik Schutte Internship 2016/Results/Trans/GENCODE/t0/gsTcell_t0_trans-eQTLs_0.05.RData")
t.0 <- trans[1:5,]
t.0 <- cbind(t.0, data.frame(t.interval=c("_t0"),origin=c("trans")))
load("~/Dropbox/Erik Schutte Internship 2016/Results/Trans/GENCODE/t10/gsTcell_t10_trans-eQTLs_0.05.RData")
t.10 <- trans[1:5,]
t.10 <- cbind(t.10, data.frame(t.interval=c("_t10"),origin=c("trans")))
load("~/Dropbox/Erik Schutte Internship 2016/Results/Trans/GENCODE/t30/gsTcell_t30_trans-eQTLs_0.05.RData")
t.30 <- trans[1:5,]
t.30 <- cbind(t.30, data.frame(t.interval=c("_t30"),origin=c("trans")))
load("~/Dropbox/Erik Schutte Internship 2016/Results/Trans/GENCODE/t180/gsTcell_t180_trans-eQTLs_0.05.RData")
t.180 <- trans[1:5,]
t.180 <- cbind(t.180, data.frame(t.interval=c("_t180"),origin=c("trans")))

### Functions
## Creates boxplots for genotypes.
# Per eQTL create a boxplot over the genotypes.
eQTLs.df.list <- list(c.0=c.0, c.10=c.10, c.30=c.30, c.180=c.180,
                      t.0=t.0, t.10=t.10, t.30=t.30, t.180=t.180)

source("~/Dropbox/Erik Schutte Internship 2016/Code/QTL-Mapping/mapping/source/prepare_df_eqtls.R")
eQTLs.df.list <- lapply(eQTLs.df.list, prepare_df_eqtls)

colnames(gene.expr) <- gsub("-","_",colnames(gene.expr))

times <- list(t0=grep("_t0",colnames(gene.expr)),
              t10=grep("_t10",colnames(gene.expr)),
              t30=grep("_t30",colnames(gene.expr)),
              t180=grep("_t180",colnames(gene.expr)))

require(reshape2)
require(ggplot2)
conversed.genotypes <- read.table("~/Dropbox/Erik Schutte Internship 2016/Data/filtered-data/snps.t.tsv",header=T)
ll <- lapply(eQTLs.df.list, function(eqtl.df) {
  apply(eqtl.df, 1, function(eqtl) {
    # Set each row as a dataframe instead of a vector
    eqtl <- data.frame(t(eqtl))
    # Position current gene in the expression file.
    gene.expression.per.interval <- gene.expr[ eqtl$gene, grep(eqtl$t.interval,colnames(gene.expr))]
    
    # Add a column to the gene expression dataframe and fill it with timepoints for later use.
    df.ge <- data.frame(gene.expression.per.interval)
    df.ge <- cbind.data.frame(df.ge, times=eqtl$t.interval)

    df.ge.melt <- melt(df.ge)
    df.ge.melt$variable <- rownames(df.ge)
    
    # Run an index on which snp corrosponds to which genotypes column.
    index = sapply(eqtl$snps, function(snp) {
      which(snp == colnames(snps))
    })
    
    # Retrieve all genotypes from all samples for the specified snps.
    genotypes <- snps[1:21,index] # because of missing sample now 1:21
    
    # Transpose and to data frame.
    genotypes <- t(data.frame(genotypes))
    
    # Melt the transposed genotypes.
    genotypes.melt <- melt(t(genotypes))
    
    # Save sample level of genotypes.
    sample.levels <- levels(genotypes.melt[,1])
    
    # Re-order genotypes on factor level.
    genotypes.melt <- genotypes.melt[order(genotypes.melt[,1]),]
    
    # Substite the sample names.
    df.ge.melt[,2] <- sapply( df.ge.melt[,2], gsub, pattern="batch[0-9]+_T", replacement = "T" )
    df.ge.melt[,2] <- sapply( df.ge.melt[,2], gsub, pattern="_t[0-9]+", replacement= "")
    
    # Bind the dataframes together.
    df.melt <- cbind.data.frame(timepoints=df.ge.melt[,1],
                                expression=df.ge.melt[,3],
                                sample=df.ge.melt[,2],
                                genotypes=genotypes.melt)
    
    # Create an order for the genotypes.
    print(as.character(eqtl$minor_major))
    major.allele = strsplit(as.character(eqtl$minor_major),"_")[[1]][2]
    print("Major allele")
    print(major.allele)
    print(strsplit(levels(df.melt$genotypes.value),"")[[1]][1])
    print(strsplit(levels(df.melt$genotypes.value),"")[[1]][2])
    if (strsplit(levels(df.melt$genotypes.value),"")[[1]][1] != major.allele |
        strsplit(levels(df.melt$genotypes.value),"")[[1]][2] != major.allele) {
      f.levels <- levels(df.melt$genotypes.value)
      f.mj <- f.levels[3]
      f.mi <- f.levels[1]
      f.he <- f.levels[2]
      if ( !is.na(f.mj) ) {
        genotype.order <- c(f.mj, f.he, f.mi)
      } else {
        genotype.order <- c(f.he, f.mi)
      }
      
      print("New levels")
      print(genotype.order)
    } else {
      f.levels <- levels(df.melt$genotypes.value)
      f.mj <- f.levels[1]
      f.he <- f.levels[2]
      f.mi <- f.levels[3]
      if ( !is.na(f.mi) ) {
        genotype.order <- c(f.mj, f.he, f.mi)
      } else {
        genotype.order <- c(f.mj, f.he)
      }
      
    }
    
    # Order the factor levels for the genotypes.
    df.melt$genotypes.value <- factor(df.melt$genotypes.value, levels=genotype.order)
    
    # Order data on factor levels.
    df.melt <- df.melt[order(df.melt$timepoints, df.melt$genotypes.value),]
    
    # Save ggplot in variable and plot.
    mi <- min( df.melt$expression ) - 1
    ma <- max(df.melt$expression ) + 1
    p <- ggplot(data=df.melt, aes(x=genotypes.value, y=expression, group=genotypes.value) ) +
      geom_boxplot(aes( fill=genotypes.value), outlier.colour = "darkgrey" ) +
      geom_point( position = position_jitter(width=0.35), colour = "darkgrey") +
      coord_cartesian( ylim = c( mi,ma ) ) +
      ggtitle( paste(eqtl$genenames, " - ", eqtl$snps, sep = ""),
               subtitle = paste( "P-value: ", eqtl$pvalue ) ) +
      theme( plot.title = element_text( size = rel(1.6), hjust = 0.5 ),
             plot.subtitle = element_text(size = rel(1), hjust = 0.5 ) ) +
      xlab(paste("Genotypes",sep="")) + ylab("Gene expression") +
      scale_fill_discrete( name="Genotypes",
                           labels=paste( names( table( df.melt$genotypes.value ) ),"(", table( df.melt$genotypes.value ), ")", sep ="") )
    p + facet_wrap( ~ timepoints, scales="free")
    
    ggsave(filename=paste("~/Dropbox/Erik Schutte Internship 2016/Images/GENCODE/",sapply(as.character(eqtl$t.interval), strsplit, split="_")[[1]][2]
              ,"/",eqtl$genenames,"_",eqtl$snps,"_",eqtl$minor_major,".pdf",sep=""),
           plot=last_plot(), device = "pdf")
    
    expr <- data.frame(gene.expr[eqtl$gene,grep(eqtl$t.interval,colnames(gene.expr))])
    geno <- data.frame(conversed.genotypes[eqtl$snps,])
    
    
    pdf(paste("~/Dropbox/Erik Schutte Internship 2016/Images/GENCODE/",sapply(as.character(eqtl$t.interval), strsplit, split="_")[[1]][2]
              ,"/",eqtl$genenames,"_",eqtl$snps,"_","correlation.pdf",sep=""), width=7, height=5)

    plot( t(geno), expr[,1], ylim=c(mi,ma), main = paste(eqtl$genenames, "-", eqtl$snps, "-", eqtl$pvalue),
          xlab = "Genotype", ylab= "Gene expression", sub = eqtl$t.interval)
    abline(lm(expr[,1] ~ t(geno) ), col="red" )
    dev.off()
    p.value <- cor.test(t(geno), expr[grep(paste(sapply(as.character(eqtl$t.interval), strsplit, split="_")[[1]][2]),rownames(expr)),] )$p.value
  })
})
