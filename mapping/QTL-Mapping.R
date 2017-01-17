############################################################
# Author: E. Schutte                                       #
# This script can map cis- trans- eQTLs and splice QTLs.   #
############################################################

# Let the user know the packages are loaded
cat("Loading packages, if not present, they are installed.\n\nThis may take a while.\n\n")

# Path for local source files.
path = file.path(getwd(),"source/")

# File names.
file.names <- dir(path, "*.R")

# Loop through file names.
for ( name in file.names ) {
  source( file.path( path, name ) )
}

## Main.
# Main function, boots script.
main <- function() {
  cat("\n\n")
  ### Call functions
  ## Settings.
  # Loads defualt settings.
  settings()
  
  ## parseOpts.
  # Parses command line options.
  parseOpts()
  
  ### Output
  ## Show output on terminal.
  # Let the user know the program started.
  cat("Starting eQTL-mapping.R..\n\nLoading files..\n\n")
  
  # Print if the files are loaded.
  cat("Files loaded:\n",as.character(args[2][[1]]),"\n",as.character(args[3]),
      "\n\nMapping Type:\n",as.character(args[1][[1]]),"\n")
  
  ### Running the directory preperation
  ## Loading files
  # Checking existing paths.
  cat("\nCheck if directorys exist..\n\n")
  
  ## Checking directories.
  # Main directory, should be universal on every system.
  global("mainDir", "~/Documents/")
  
  # Sub directory.
  global("subDir", "R/")
  
  ### Prepare data
  ## prepare data files for mapping.
  # Call data_prep to prepare data.
  data_prep()
  
  ## Check Type.
  # Check the mapping type.
  if ( MAPPING_TYPE == 'basic') {
    
    # Create directory for the Basic results.
    createDir("Basic")
    
    # Map eQTLs.
    eQTL()
    
  } else if ( MAPPING_TYPE == 'ct' ) {
    
    # Create directory for the Cis results.
    createDir("Cis")
    
    # Create directory for the Trans results.
    createDir("Trans")
    
    # Map cis- trans-eQTLs.
    ctQTL()
    
  } else if ( MAPPING_TYPE == 's' ) {
    
    # Create directory for the Splice results.
    createDir("Splice")
    
    # Map splice QTLs.
    sQTL()
    
  } else {
    error("No mapping type accepted", "Something unexpected")
  }
}

main()

# ## CRAN
# # Check if all libraries are installed, if not install them.
# packages <- c("MatrixEQTL","devtools","plyr","gtools","dplyr")
# if ( length( setdiff( packages, rownames( installed.packages() ) ) ) > 0) {
#   suppressMessages(
#     install.packages( setdiff( packages, rownames( installed.packages() ) ), repos = "http://cran.xl-mirror.nl" )  
#   )
# }
# 
# ## BiocLite
# # Load BiocLite packages.
# bioclite.packages <- c("vegan", "Rsamtools", "qvalue","Biostrings")
# if( length( setdiff( bioclite.packages, rownames( installed.packages() ) ) > 0 ) ) {
#   suppressMessages(
#     biocLite( setdiff( bioclite.packages, rownames( installed.packages() ) ) )
#   )
# }
# 
# ## Libraries
# # Load Matrix eQTL library.
# library("MatrixEQTL")
# library("devtools")
# if ( length( dev_packages() ) < 1 ) {
#   devtools::install_github("guigolab/sQTLseekeR")
# } 
# library("sQTLseekeR")
# library("gtools")
# library("plyr")
# library("dplyr")