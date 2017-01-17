## Settings.
# Set all needed settings.
settings <- function() {
  
  ## Command line arguments.
  # Clear the current environment.
  rm(list=ls())
  
  # Makes command line arguments availble for this script.
  global("args",commandArgs(trailingOnly = T))
  
  # Command line Test 1:
  #global("args",list("-h"))
  
  # Command line Test 2:
  #global("args",list("basic","~/Dropbox/Erik/clones_CeD_genotypes_noDuplicates.Rdata","~/Dropbox/Erik/count_all_batch_times_sample.names_96Samples_noCD8_leiden.Rdata"))
  
  # Command line Test 3:
  #global("args", list("basic","~/Dropbox/Erik/clones_CeD_genotypes_noDuplicates.Rdata","~/Dropbox/Erik/count_all_batch_times_sample.names_96Samples_noCD8_leiden.Rdata","-v"))
  
  # Command line Test 4:
  #global("args", list("basic","~/Dropbox/Erik/count_all_batch_times_sample.names_96Samples_noCD8_leiden.Rdata","-v"))
  
  # Command line Test 5:
  #global("args",list("basic","~/Dropbox/Erik/count_all_batch_times_sample.names_96Samples_noCD8_leiden.Rdata"))
  
  # Command line Test 6:
  #global("args",list("~/Dropbox/Erik/clones_CeD_genotypes_noDuplicates.Rdata","~/Dropbox/Erik/count_all_batch_times_sample.names_96Samples_noCD8_leiden.Rdata","-v"))
  
  # Command line Test 7:
  #global("args", list("~/Dropbox/Erik/clones_CeD_genotypes_noDuplicates.Rdata","~/Dropbox/Erik/count_all_batch_times_sample.names_96Samples_noCD8_leiden.Rdata"))
  
  # Command line Test 8:
  # global("args",list("s","~/Dropbox/Erik/clones_CeD_genotypes_noDuplicates.Rdata","~/Dropbox/Erik/expr_vst_condition_patient_rmBatch_88samples.Rdata"))
  
  #load("~/Dropbox/Erik/clones_CeD_genotypes_noDuplicates.Rdata")
  #load("~/Dropbox/Erik/count_all_batch_times_sample.names_96Samples_noCD8_leiden.Rdata")
  #load("~/Dropbox/Erik/expr_vst_condition_patient_rmBatch_88samples.Rdata")
  
  # Set mapping flag.
  global("MAPPING_TYPE", "")
  global("MAPPING_OPTS", c("basic","ct","s"))
  
  # Set verbose flag.
  global("VERBOSE", FALSE)
}
