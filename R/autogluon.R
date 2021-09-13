#' Function to predict the GEP-subgroups using Autogluon model
#'
#' @param x Normalized gene expression matrix or its file path
#' @param exp_type deseq2 or log2tpm
#' @export
#' @examples
#' x <- system.file("extdata", "tcga_log2tpm_test.csv", package = "clinaml")
#' res <- autogluon_pred(x, exp_type = "log2tpm")
#' x2 <- system.file("extdata", "tcga_deseq2_test.csv", package = "clinaml")
#' res2 <- autogluon_pred(x2, exp_type = "deseq2")
autogluon_pred <- function (x = NULL, exp_type = "log2tpm") {
  if (is.character(x) && file.exists(x)) {
    x <- fread(x, data.table = FALSE, fill = TRUE)
  } else if (is.character(x)) {
    return(NULL)
  }
  x <- check_gep_dat(x)
  model_dat <- get_model_data()[["Autogluon"]]
  model_dat <- model_dat[str_detect(names(model_dat), exp_type)]
  script <- system.file("extdata", "autogluon.py", package = "clinaml")
  tmp_in <- tempfile(fileext = ".csv")
  fwrite(as.data.frame(x), tmp_in, sep=",", row.names = FALSE, quote = FALSE)
  res_list <- list(pred = NULL, summary = NULL)
  out <- list()

  condaenv <- options("clinaml.condaenv")[[1]]
  condaenv <- ifelse(is.null(condaenv), "clinaml", condaenv)
  use_condaenv(condaenv, required = TRUE)
 
  for (i in names(model_dat)) {
    tmp_out <- tempfile(fileext = ".csv")
    cmd <- sprintf("python %s %s %s %s", script, model_dat[i], tmp_in, tmp_out)
    cat(cmd, sep = "\n")
    system(cmd)
    out[[i]] <- tmp_out
    res_list$pred[[i]][[1]] <- fread(tmp_out, data.table = FALSE)[, -1]
    pred <- res_list$pred[[i]][[1]]
    res_list$pred[[i]][[2]] <- colnames(pred)[apply(pred, 1, function(x) {which(x == max(x))})]
  }
  res_list <- get_pred_summary(res_list, names(model_dat))
}
