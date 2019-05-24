.onLoad <- function(libname, pkgname) {
  # make data set names global to avoid CHECK notes
  utils::globalVariables(c("True", "Predicted", "CLUSTER", "inside", "Cell_type",
                           "Cluster", "MEAN", "MEDIAN", "MAX"))
  invisible()
}
