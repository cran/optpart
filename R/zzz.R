.First.lib <- function(lib, pkg) {
  library.dynam("optpart", pkg, lib)
  require(MASS)
  require(stats)
  require(cluster)
  require(plotrix)
}
