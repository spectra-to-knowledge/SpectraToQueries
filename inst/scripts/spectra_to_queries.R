start <- Sys.time()

message("This program extracts diagnostic ions from groups of spectra.")
message("Authors: Adriano Rutz, Sarah Szwarc")

SpectraToQueries::spectra_to_queries()

end <- Sys.time()

message("Script finished in ", format(end - start))
