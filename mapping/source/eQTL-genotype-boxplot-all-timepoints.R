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

### ALL TIMEPOINTS
## Load data files for eqtl plots
# Cis eqtls NOTE: interactive or just all timepoints.
#load("~/Dropbox/Erik Schutte Internship 2016/Results/Cis/GENCODE/gsTcell_all_timepoints_cis-eQTLs_0.05.RData")
load("~/Dropbox/Erik Schutte Internship 2016/Results/Cis/GENCODE/gsTcell_all_timepoints_interactive-cis-eQTLs_0.05.RData")
c.all <- cis[1:10,]
c.all <- cbind(c.all, data.frame(t.interval=c("all"),origin=c("interactive_cis")))

# Trans eqtls NOTE: interactive or just all timepoints.
#load("~/Dropbox/Erik Schutte Internship 2016/Results/Trans/GENCODE/gsTcell_all_timepoints_trans-eQTLs_0.05.RData")
load("~/Dropbox/Erik Schutte Internship 2016/Results/Trans/GENCODE/gsTcell_all_timepoints_interactive-trans-eQTLs_0.05.RData")
t.all <- trans[1:10,]
t.all <- cbind(t.all, data.frame(t.interval=c("all"),origin=c("interactive_trans")))


### Functions
## Creates boxplots for genotypes.
# Per eQTL create a boxplot over the genotypes.
eQTLs.df.list <- list(c.all=c.all, t.all=t.all)

source("~/Dropbox/Erik Schutte Internship 2016/Code/QTL-mapping/LncRNA_Explorer_Scripts/prepare_df_eqtls.R")
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
    gene.expression.all.samples <- gene.expr[ which( eqtl$gene == rownames(gene.expr) ), ]

    # Add a column to the gene expression dataframe and fill it with timepoints for later use.
    df.ge <- data.frame(gene.expression.all.samples)
    df.ge <- cbind.data.frame(df.ge, times=0)
    df.ge[times$t0,2] <- "t0"
    df.ge[times$t10,2] <- "t10"
    df.ge[times$t30,2] <- "t30"
    df.ge[times$t180,2] <- "t180"
    df.ge.melt <- melt(df.ge)
    df.ge.melt[,2] <- rownames(df.ge)

    # Run an index on which snp corrosponds to which genotypes column.
    index = sapply(eqtl$snps, function(snp) {
      which(snp == colnames(snps))
    })

    # Retrieve all genotypes from all samples for the specified snps.
    genotypes <- snps[1:21,index] # because of missing sample now 1:21

    # Transpose and to data frame.
    genotypes <- t(data.frame(genotypes))

    # Multiply the genotypes 4 times.
    genotypes <- cbind.data.frame(genotypes, genotypes, row.names = NULL) # 44x / 42x
    genotypes <- cbind.data.frame(genotypes, genotypes, row.names = NULL) # 88x / 84x

    # Melt the transposed genotypes.
    genotypes.melt <- melt(t(genotypes))

    # Save sample level of genotypes.
    sample.levels <- levels(genotypes.melt[,1])

    # Re-order genotypes on factor level.
    genotypes.melt <- genotypes.melt[order(genotypes.melt[,1]),]

    # Substite the sample names.
    df.ge.melt[,2] <- sub("batch[0-9]+__T","T",df.ge.melt[,2])
    df.ge.melt[,2] <- sub("__t[0-9]+","",df.ge.melt[,2])

    # Bind the dataframes together.
    df.melt <- cbind.data.frame(timepoints=df.ge.melt[,1],
                                expression=df.ge.melt[,3],
                                sample=df.ge.melt[,2],
                                genotypes=genotypes.melt)

    # Create an order for the timepoints.
    timepoints.order <- c("t0","t10","t30","t180")

    # Order the factor levels for the timepoints.
    df.melt$timepoints <- factor(df.melt[,1], levels=timepoints.order)
    print(eqtl)
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
                          labels=paste( names( table( df.melt$genotypes.value ) ),"(", table( df.melt$genotypes.value )/4, ")", sep ="") )
    p + facet_wrap( ~ timepoints, scales="free")

    ggsave(filename=paste("~/Dropbox/Erik Schutte Internship 2016/Images/GENCODE/interactive/",eqtl$genenames,"_",eqtl$snps,"_",eqtl$minor_major,"_",eqtl$origin,".pdf",sep=""),
           plot=last_plot(), device = "pdf")

    expr <- data.frame(gene.expr[which(eqtl$gene == rownames(gene.expr)),])
    geno <- data.frame(conversed.genotypes[which(rownames(conversed.genotypes) == eqtl$snps),])


    pdf(paste("~/Dropbox/Erik Schutte Internship 2016/Images/GENCODE/interactive/",eqtl$genenames,"_",eqtl$snps,"_",eqtl$minor_major,"_",eqtl$origin,"_","correlation.pdf",sep=""),
        width=7, height=5)
    par(mfrow=c(2,2))
    count = 1
    for ( time in times ) {
      plot( t(geno), expr[time,], ylim=c(mi,ma) , main = paste("Gene: ", eqtl$genenames, "\nSNP: ", eqtl$snps ),
            xlab = "Genotype", ylab = "Gene expression", sub = names(times)[count]  )
      abline(lm(expr[time,] ~ t(geno) ), col="red" )
      count = count + 1
    }
    dev.off()
    p.value <- cor.test(t(geno), expr[grep("t0",rownames(expr)),])$p.value

    # p.other <- as.double(as.character(eqtl$pvalue))
    # print("Pvalue and class cor:")
    # print(p.value)
    # print(class(p.value))
    # print("Pvalue and class eqtl:")
    # print(p.other)
    # print(class(p.other))
    # print("Are they the same?")
    # print(p.value == p.other)
    # print("eqtl pvalue as factor")
    # print(eqtl$pvalue)
    # imp.genotypes.snp <- conversed.genotypes[which(eqtl$snps == rownames(conversed.genotypes)),]
    # ge.snp <- gene.expr[ which( eqtl$gene == rownames(gene.expr) ), ]
    # ge.snp <- t(matrix(ge.snp))
    # imp.genotypes.snp <- cbind(imp.genotypes.snp, imp.genotypes.snp) # 44x
    # imp.genotypes.snp <- cbind(imp.genotypes.snp, imp.genotypes.snp) # 88x
    # p.value <- cor.test(as.matrix(imp.genotypes.snp), ge.snp, method="spearman")$p.value
    # if (p.value != as.double(eqtl$pvalue)) {
    #   cat("Cor.test: ",p.value, " not equal to MatrixEQTL: ", as.double(eqtl$pvalue))
    # } else {
    #   cat("P-values match!")
    # }
  })
})


