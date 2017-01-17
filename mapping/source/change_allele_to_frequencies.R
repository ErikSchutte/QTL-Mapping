## Change_allele_to_frequencies.
# Function to identify the genotype and determine wheter it is major, minor or heterozygote.
change_allele_to_frequencies <- function (x,allele.list) {
  if ( identical(x,allele.list[[1]]$genotype) ) {
    # The genotype is major, return 0.
    return(0)
    
  } else if ( identical(allele.list[[4]]$genotype,x) ) {
    # The genotype is NA or 00, return NA.
    return('NA')
    
  } else if ( identical(x,allele.list[[3]]$genotype[1]) | identical(x,allele.list[[3]]$genotype[2]) ) {
    # If the genotype is heterozygote, return 1.
    return(1)
    
  } else {
    # If the genotype is one of the minor alleles, return 2.
    return(2)
  }
}