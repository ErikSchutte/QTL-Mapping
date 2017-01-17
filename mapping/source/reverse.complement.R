## reverse.complement.
# Function that takes the reverse complement of the Major allele in case there are only major allels and hetrozygous alles,
# we need to know what the minor allele is for the database.
reverse.complement <- function (major.allele) {
  require(Biostrings)
  major.allele <- DNAString(major.allele[1])
  major.allele <- reverseComplement(major.allele)
  major.allele <- strsplit(as.character(major.allele)," ")
  return(major.allele[[1]][1])
}