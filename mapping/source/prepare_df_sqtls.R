prepare_df_sqtls <- function(sqtl) {
  tmp.df <- NULL
  minor_allele <- c()
  minor_alleles <- sapply(sqtl$SNP, function(snp) {
    alleles <- alleles[which(alleles[,2]==snp),1]
    return(alleles)
  })
  tmp.df <- cbind(sqtl)
  tmp.df <- cbind.data.frame(tmp.df, data.frame(genotype=minor_alleles))
  tmp.df <- tmp.df[c("gene","SNP","t.stat","p.value","FDR","beta","genotype")] # Implement genenames when added.
  colnames(tmp.df) <- c("gene", "snps", "statistic", "pvalue", "FDR", "beta", "minor_major") # Implement genenames when added.

  return(tmp.df)
}