#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(format(utils::citation(pkgname)))
}
