## parseOpts.
# Parse command line options and test them.
parseOpts <- function() {
  
  # If length of arguments is 1.
  if ( length(args) == 1 ) {
    
    # Check if the first argument is a help flag.
    if ( length(args[1][[1]]) != 0 && args[1][[1]] == '-h' | length(args[1][[1]]) != 0 && args[1][[1]] == '--help')  {
      help()
      
    } 
    # Check if the first argument is not empty or the help flag. I do this twice because I sometimes see that an argument that is not -h
    # it sometimes gets accepted.
    else if ( !is.null(args[1][[1]]) && args[1][[1]] == '-h' | args[1][[1]] == '--help' ) {
      help()
    }
    # if it is the s mapping option set it.
    else if ( args[1][[1]] %in% MAPPING_OPTS && args[1][[1]] == 's') {
      global( "MAPPING_TYPE", as.character(args[1][[1]]) )
    }
    # If it is not a help flag, throw error unknown flag.
    else {
      error( "Unknown flag", args[1][[1]] )
    }
    
    # If length of arguments is 2 or 3.
  } else if ( length(args) == 2 | length (args) == 3 ) {
    
    # If the length of the first argument is not 0 or the first argument is not null.
    if ( length(args[1][[1]]) != 0 | !is.null(args[1][[1]]) ) {
      
      # If the first argument is in MAPPING_OPTS.
      if ( args[1][[1]] %in% MAPPING_OPTS ) {
        
        # Set first argument as global variable MAPPTING_TYPE.
        global( "MAPPING_TYPE", as.character(args[1][[1]]) )
        
        # If the first argument is not in MAPPING_OPTS throw error invalid mapping type.
      } else {
        error( "Invalid mapping type.",args[1][[1]] )
      }
    }
    
    # If argument 2 is zero or null and argument 3 is zero or null.
    if ( length(args[2][[1]]) != 0  && length(args[3][[1]]) != 0  ) {
      #| !is.null(args[2][[1]])  | !is.null(args[3][[1]])
      
      # Check if file argument 2 and argument 3 exist.
      if ( file.exists( file.path(args[2][[1]]) ) && file.exists( file.path(args[3][[1]]) ) ) {
        
        # Create temp space.
        temp.space <- new.env()
        
        # Load data into temp space.
        snps.temp <- load(args[2][[1]], temp.space)
        GE.temp <- load(args[3][[1]], temp.space)
        
        # Save data to current environment.
        snps <- get(snps.temp, temp.space) # Genotype file.
        GE <- get(GE.temp, temp.space) # Gene expression file.
        
        # Save data globally.
        global("GE", GE)
        global("snps", snps)
        
        # Store all the names in the global environment.
        #nms <- ls(pattern = "", envir = temp.space)
        #L <- sapply(nms, get, envir = .GlobalEnv, simplify = FALSE)
        
        # Sort out the matrixes.
        #L <- as.matrix(L[which(sapply(L,is.matrix) == T)])
        
        # If files don't exist, throw error.
      } else {
        
        # If the -v, --verbose flag is used as file input, throw an error.
        if ( args[2][[1]] == '-v' | args[2][[1]] == '--verbose' | args[3][[1]] == '-v' | args[3][[1]] == '--verbose' ) {
          error("-v, --verbose flag used as file input",args[2:3])
        }
        
        error("File/files does/do not exist",args[2:3])
      }
      
      # If argument 2 and 3 are empty, throw error to provide files.
    } else {
      error("Provide two input files",args[2:3])
    }
    
    # If length of args is 4.  
  } else if ( length (args) == 4 ) {
    
    # If the length of the first argument is not null or the first argument is not null.
    if ( length(args[1][[1]]) != 0 | !is.null(args[1][[1]]) ) {
      
      # If the first argument is in MAPPING_OPTS.
      if (args[1][[1]] %in% MAPPING_OPTS) {
        
        # Set the first argument as the MAPPING_TYPE.
        global("MAPPING_TYPE", as.character(args[1][[1]]))
        
        # If the first argument is not in MAPPING_OPTS, throw error invalid mapping type.  
      } else {
        error("Invalid mapping type.",args[1][[1]])
      }
    }
    
    # If the length of the 2nd argument is not null or the 2nd argument is not null,
    # and if the length of the 3th argument is not null or the 3th argument is not null.
    if ( length(args[2][[1]]) != 0  && length(args[3][[1]]) != 0  ) {
      
      # If the file as 2nd argument exists and the file as 3th argument exists.
      if( file.exists( file.path(args[2][[1]]) ) && file.exists( file.path(args[3][[1]]) ) ) {
        
        # Create temp space.
        temp.space <- new.env()
        
        # Load data into temp space.
        snps.temp <- load(args[2][[1]], temp.space)
        GE.temp <- load(args[3][[1]], temp.space)
        
        # Save data to current environment.
        snps <- get(snps.temp, temp.space) # Genotype file.
        GE <- get(GE.temp, temp.space) # Gene expression file.
        
        # Save data globally.
        global("GE", GE)
        global("snps", snps)
        
        # Store all the names in the global environment.
        #nms <- ls(pattern = "", envir = temp.space)
        #L <- sapply(nms, get, envir = .GlobalEnv, simplify = FALSE)
        
        # Sort out the matrixes.
        #L <- as.matrix(L[which(sapply(L,is.matrix) == T)])
        
        # If files don't exist, throw error.
      } else {
        
        # If the -v, --verbose flag is used as file input, throw an error.
        if ( args[2][[1]] == '-v' | args[2][[1]] == '--verbose' | args[3][[1]] == '-v' | args[3][[1]] == '--verbose' ) {
          error("-v, --verbose flag used as file input",args[2:3])
        }
        error("File/files does/do not exist",args[2:3])
      }
      
      # If the length of the 2nd argument is null or the 2nd argument is null,
      # and if the length of the 3th argument is null or the 3th argument is null,
      # throw an error to provide two input files.
    } else {
      error("Provide two input files",args[2:3])
    }
    
    # If the length of the 4th argument is not null or the 4th argument is not null.
    if ( length(args[4][[1]]) != 0 | !is.null(args[4][[1]]) ) {
      
      # If the 4th argument equals -v or --verbose.
      if ( args[4][[1]] == "-v" | args[4][[1]] == "--verbose" ) {
        
        # Set the VERBOSE variable to TRUE.
        global("VERBOSE", TRUE)
        
        # Throw an error of an unkown flag.  
      } else {
        error("Unknown flag",args[4][[1]])
      }
    }
    
    # If no arguments are provided throw an error.
  } else {
    error("No arguments provided",args)
  }
}