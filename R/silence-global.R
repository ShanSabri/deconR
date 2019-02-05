.onLoad <- function(libname, pkgname) {
  # make data set names global to avoid CHECK notes
  utils::globalVariables("True")
  utils::globalVariables("Predicted")
  utils::globalVariables("CLUSTER")
  utils::globalVariables("inside")
  utils::globalVariables("Cell_type")

  invisible()
}
