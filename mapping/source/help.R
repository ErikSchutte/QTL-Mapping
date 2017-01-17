## Help.
# display help message.
help <- function() {
  cat("
      QTL-mapping Script
      
      This script maps genotypes from samples to expression values, locating basic-, cis-, trans-
      and s-QTLs. Important to note; The dimensions from the genotype file and the gene expression
      file have to be from the same length, or a multiplication of that same length (given time
      points).
      
      !!The data has to be saved in an *.Rdata image, that is then loaded into the script.!!
      
      Arguments:
      --help, -h                   - Displays this help message.
      
      map_type                     - The sort of eQTL Analysis you want performed
      |- Options:
      |- 'basic', basic eQTL analysis.
      |- 'ct', cis- trans-eQTL analysis.
      |- 's', sQTL analysis (NOT IMPLMENTED).
      
      genotype_file                - A File containing the genotypes for merged alleles per sample (column),
      per gene (row).
      
      gene_expression_file         - A File containing the gene expressino data per sample (column)
      per gene (row).
      
      --verbose, -v                - A Verbose option. Shows output of steps taken. You might want to save
      the std.out to a file ( see example ).
      
      Examples:
      - Texual;
      |-Rscript --vanilla eQTL-mapping.R <map_type> <genotype_file> <gene_expression_file> [--verbose, -v]
      
      - Script call;
      |-Rscript --vanilla eQTL-mapping.R basic myGenotypeFile.Rdata myGeneExpressionFile.Rdata
      |
      |- With verbose output;
      |  |-Rscript --vanilla eQTL-mapping.R basic myGenotypeFile.Rdata myGeneExpressionFile.Rdata -v
      |
      |- With verbose output, saved to a file;
      |-Rscript --vanilla eQTL-mapping.R basic myGenotypeFile.Rdata myGeneExpressionFile.Rdata -v > output.txt
      
      ")
  q(save="no")
}