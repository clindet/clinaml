#' Function to predict the GEP-subgroups using XGBoost model
#'
#' @param x Normalized gene expression matrix or its file path
#' @param exp_type deseq2 or log2tpm
#' @param model_counter counter of reads for modeling
#' @param model_sva sva processed or not for modeling
#' @param model_folds model fold
#' @export
#' @examples
#' x <- system.file("extdata", "tcga_log2tpm_test.csv", package = "clinaml")
#' res <- xgboost_pred(x, exp_type = "log2tpm")
#' x2 <- system.file("extdata", "tcga_deseq2_test.csv", package = "clinaml")
#' res2 <- xgboost_pred(x2, exp_type = "deseq2")
xgboost_pred <- function (x = NULL, exp_type = "log2tpm",
  model_counter = c("featurecounts", "htseq", "kallisto", "salmon"),
  model_sva = c("yes", "no"), model_folds = 1:10) {
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

  exp_type <- ifelse(exp_type == "log2tpm", "log2tpm", "deseq2vsd") 
  sva_prefix <- c(yes="after_sva", "no"="before_sva")
  sva_prefix <- sva_prefix[model_sva]
  x1 <- paste(model_counter, exp_type, sep = "_")
  x2 <- unlist(sapply(x1, function(x) {paste0(x, "_", sva_prefix, "_fold")}))
  mnames <- unlist(sapply(x2, function(x) paste0(x, model_folds)))
  mnames <- as.character(mnames)
  print(mnames)
  for (i in mnames) {
    xgb <- readRDS(model_dat[[i]])
    res_list$pred[[i]] <- prob_mat_to_label(xgb, xgb_test, nrow(x), 6)
  }
  res_list <- get_pred_summary(res_list, names(model_dat))
  return(res_list)
}
