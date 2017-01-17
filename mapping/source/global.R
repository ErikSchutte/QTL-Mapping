## Global settings.
# Stores settings globally.
global <- function(x,y) {
  assign(x,y,envir=.GlobalEnv)
}