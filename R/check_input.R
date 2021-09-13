#' Function to check whether the input fits all features for
#' GEP-defined subgroups prediction
#'
#' @param x Normalized gene expression matrix,
#' the 1st column is gene id, ensemble id, or gene symbol
#' @export
#' @examples
#' check_gep_dat(NULL)
#' @return logical
check_gep_dat <- function (x) {
  colnames(x) <- str_replace_all(colnames(x), fixed("."), "-")
  idx <- c()
  for (i in 1:ncol(x)) {
    idx <- c(idx, is.character(x[,i]))
  }
  return(x[,!idx])
}
