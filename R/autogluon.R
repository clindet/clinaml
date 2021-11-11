#' Function to predict the GEP-subgroups using Autogluon model
#'
#' @param x Normalized gene expression matrix or its file path
#' @param exp_type deseq2 or log2tpm
#' @param model_counter counter of reads for modeling
#' @param model_sva sva processed or not for modeling
#' @param model_folds model fold
#' @export
#' @examples
#' x <- system.file("extdata", "tcga_log2tpm_test.csv", package = "clinaml")
#' res <- autogluon_pred(x, exp_type = "log2tpm")
#' x2 <- system.file("extdata", "tcga_deseq2_test.csv", package = "clinaml")
#' res2 <- autogluon_pred(x2, exp_type = "deseq2")
autogluon_pred <- function (x = NULL, exp_type = "log2tpm", 
  model_counter = c("featurecounts", "htseq", "kallisto", "salmon"),
  model_sva = c("yes", "no"), model_folds = 1:10) {
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

  #condaenv <- options("clinaml.condaenv")[[1]]
  #condaenv <- ifelse(is.null(condaenv), "clinaml", condaenv)
  #use_condaenv(condaenv)
  exp_type <- ifelse(exp_type == "log2tpm", "log2tpm", "deseq2vsd") 
  sva_prefix <- c(yes="after_sva", "no"="before_sva")
  sva_prefix <- sva_prefix[model_sva]
  x1 <- paste(model_counter, exp_type, sep = "_")
  x2 <- unlist(sapply(x1, function(x) {paste0(x, "_", sva_prefix, "_fold")}))
  mnames <- unlist(sapply(x2, function(x) paste0(x, model_folds)))
  for (i in as.character(mnames)) {
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
