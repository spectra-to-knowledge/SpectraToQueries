start <- Sys.time()

message("This program extract diagnostic ions from groups of spectra.")
message("Authors: \n", "AR")
message("Contributors: \n", "...")

SpectraToQueries::spectra_to_queries()

end <- Sys.time()

message("Script finished in ", format(end - start))
