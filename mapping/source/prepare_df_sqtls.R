prepare_df_sqtls <- function(sqtl) {
  # Different components from cluster. [1]=chr, [2]=intron_start, [3]=intron_end, [4]=cluster_id
  # names <- sapply(sqtl[,2], function(pos){
  #   components <- strsplit(as.character(pos),":")[[1]]
  #   chr=components[1] 
  #   chr="chr1"
  #   intron_start=components[2]
  #   intron_start=7913180
  #   intron_end=7913367
  #   intron_end=components[3]
  #   cluster=components[4]
  #   
  #   match.chr=grep(genepos[,2],pattern=chr)
  #   row = genepos[match.chr,]
  #   cat("intron start position",intron_start,"gene start",row[,3],"\nintron stop position",intron_end,"gene stop",row[,4],"\n")
  #   print(row[which(intron_start >= row[,3] & intron_end <= row[,4]),1])
  # })
  # 

  tmp.df <- NULL
  minor_allele <- c()
  minor_alleles <- sapply(as.character(sqtl[,1]), function(snp) {
    alleles <- alleles[which(alleles[,2]==snp),1]
    return(alleles)
  })
  tmp.df <- cbind(sqtl)
  tmp.df <- cbind.data.frame(tmp.df, data.frame(genotype=minor_alleles))
  tmp.df <- tmp.df[c("gene","snps","statistic","pvalue","FDR","beta","genotype")] # Implement genenames when added.
  colnames(tmp.df) <- c("gene", "snps", "statistic", "pvalue", "FDR", "beta", "minor_major") # Implement genenames when added.

  return(tmp.df)
}