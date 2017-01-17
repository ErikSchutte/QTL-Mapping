## Error handling.
# Handles error messages, exits script.
error <- function(error,input) {
  tmp <- ""
  if ( length(input) >= 2 ) {
    for ( i in 1:length(input) ) {
      tmp <- paste(as.character(input[[i]]), sep="")
    }
  }
  input <- tmp
  cat("
      Something went wrong:
      
      Your input:",
      input,"is probably invalid.
      Reason:",
      error,
      
      "
      Please check for spelling errors, and if the correct file/option is entered as input.
      ")
  q(save="no")
}