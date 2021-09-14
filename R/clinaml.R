#' The 'clinaml' package provides R functions to predict the
#' gene expression profiling (GEP)-defined subgroups and the GEP-plus risk scores.
#'
#' @author
#' Li Jianfeng \url{lee_jianfeng@sjtu.edu.cn}
#' @seealso
#' Useful links:
#'
#' \url{https://github.com/clindet/clinaml}
#'
#' Report bugs at \url{https://github.com/clindet/clinaml/issues}
#'
#' @docType package
#' @import xgboost
#' @import DESeq2
#' @import tximport
#' @import BiocParallel
#' @importFrom stringr str_detect str_replace_all str_split fixed
#' @importFrom data.table fread fwrite
#' @importFrom xgboost xgb.DMatrix
#' @importFrom utils browseURL
#' @importFrom reticulate use_condaenv
#' @importFrom parallel detectCores
#' @importFrom R.utils gunzip
#' @name clinaml
#'
NULL

.onAttach <- function(libname, pkgname) {
  op <- options()
  Sys.setenv(R_TESTS = "")
  invisible()
}
