#' Function to predict the GEP-subgroups using XGBoost model
#'
#' @param x Normalized gene expression matrix or its file path
#' @param exp_type deseq2 or log2tpm
#' @export
#' @examples
#' x <- system.file("extdata", "tcga_log2tpm_test.csv", package = "clinaml")
#' res <- xgboost_pred(x, exp_type = "log2tpm")
#' x2 <- system.file("extdata", "tcga_deseq2_test.csv", package = "clinaml")
#' res2 <- xgboost_pred(x2, exp_type = "deseq2")
xgboost_pred <- function (x = NULL, exp_type = "log2tpm") {
  if (is.character(x) && file.exists(x)) {
    x <- fread(x, data.table = FALSE, fill = TRUE)
  } else if (is.character(x)) {
    return(NULL)
  }
  x <- check_gep_dat(x)
  model_dat <- get_model_data()[["XGBoost"]]
  model_dat <- model_dat[str_detect(names(model_dat), exp_type)]
  xgb_test <- xgb.DMatrix(as.matrix(x))

  res_list <- list(pred = NULL, summary = NULL)
  for (i in names(model_dat)) {
    xgb <- readRDS(model_dat[[i]])
    res_list$pred[[i]] <- prob_mat_to_label(xgb, xgb_test, nrow(x), 6)
  }
  res_list <- get_pred_summary(res_list, names(model_dat))
  return(res_list)
}
