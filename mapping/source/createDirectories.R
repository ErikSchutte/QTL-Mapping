## createDir.
# creates direcotry for mapping types.
createDir <- function(d) {
  
  # If directory exists..
  if ( dir.exists( file.path(mainDir, subDir, d, fsep="") ) ) {
    
    # .. show the user the directory already exists.
    cat("Directory: \'", file.path(mainDir, subDir, d, fsep=""), "\' exists.\n", sep = "") 
  } 
  
  # If direcotry does not exists..
  else {
    
    # .. show the user the directory is being created.
    cat("Creating directory: \'", file.path(mainDir, subDir, d, fsep=""),"\'\n",sep="")
    
    # Create the directory.
    dir.create( file.path(mainDir, subDir, d, fsep=""), recursive = T )
  }
  cat("\n\n")
}
