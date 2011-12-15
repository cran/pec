.First.lib <- function(lib,pkg) {
  ## require(survival)
  ## require(prodlim)
  library.dynam("pec",pkg,lib)
}
.Last.lib <- function(lib)
  library.dynam.unload("pec",lib)
