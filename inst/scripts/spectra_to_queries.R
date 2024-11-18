start <- Sys.time()

pkgload::load_all()

message("This program extract diagnostic ions from groups of spectra.")
message("Authors: \n", "AR")
message("Contributors: \n", "...")

spectra_to_queries()

end <- Sys.time()

message("Script finished in ", format(end - start))
