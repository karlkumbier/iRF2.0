.onLoad <- function(libname, pkgname) {
  message('Compiling Rcpp function...')
  Rcpp::sourceCpp("hack/ancestorPath.cpp")
}
