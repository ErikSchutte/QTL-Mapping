prepare_df_eqtls <- function(eqtl) {
  tmp.df <- NULL
  # print(eqtl)
  genes <- eqtl$gene
  genenames <- c()
  for ( gene in genes ) {
    associatedgenename.index = biomart[which( gene == biomart$gene_id ),2]
    genenames <- c(genenames, associatedgenename.index)
  }
  df.snps <- eqtl$snps
  
  minor_allele <- c()
  for ( snp in df.snps ) {
    associated_alleles.index = snps[22,which( snp == colnames(snps) )]
    minor_allele <- c(minor_allele, associated_alleles.index)
  }
  
  tmp.df <- cbind(eqtl, genenames=biomart[genenames,2])
  tmp.df <- cbind.data.frame(tmp.df, data.frame(genotype=minor_allele))
  tmp.df <- tmp.df[c("gene","genenames","snps","statistic","pvalue","FDR","beta","genotype","t.interval","origin")]
  colnames(tmp.df) <- c("gene", "genenames", "snps", "statistic", "pvalue", "FDR", "beta", "minor_major","t.interval","origin")
  return(tmp.df)
}