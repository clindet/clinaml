#' Function to get the path of model data. Root path is determined by
#' option var clinaml.MODEL_DATA or "~/.clinaml/model"
#'
#' @export
#' @examples
#' get_model_data()
get_model_data <- function () {
  #options(clinaml.MODEL_DATA="D:/paper/03-todo/aml/analysis/pred_models/GEP.subgroups.model")
  #options(clinaml.MODEL_DATA="/public/home/lijf/env/model/clinaml")
  x <- options("clinaml.MODEL_DATA")[[1]]
  if (is.null(x)) {
    x <- "~/.clinaml/model"
  }
  x2 <- list(XGBoost=NULL, Autogluon=NULL)
  autogluon_model_pth <- system.file("extdata", "autogluon_model_file", package = "clinaml")
  autogluon_model <- read.csv(autogluon_model_pth, sep = ":", header = FALSE)
  autogluon_model_map <- autogluon_model$V2
  names(autogluon_model_map) <- autogluon_model$V1
  for (m in c("XGBoost", "Autogluon")) {
    for (tp in c("deseq2vsd_after_sva", "deseq2vsd_before_sva", "log2tpm_after_sva", "log2tpm_before_sva")) {
      for (counter in c("featurecounts", "htseq", "salmon", "kallisto")) {
        for (f in 1:10) {
          k <- sprintf("%s_%s_fold%s", counter, tp, f)
          if (m == "XGBoost") {
            x3 <- sprintf("%s/%s/%s.model.rds", x, m, k)
            x2[["XGBoost"]] <- c(x2[["XGBoost"]], x3)
            names(x2[["XGBoost"]])[length(x2[["XGBoost"]])] <- k
          } else {
            x3 <- sprintf("%s/%s/%s/AutogluonModels/%s", x, m, k, autogluon_model_map[k])
            x2[["Autogluon"]] <- c(x2[["Autogluon"]], x3)
            names(x2[["Autogluon"]])[length(x2[["Autogluon"]])] <- k
          }
        }
      }
    }
  }
  return(x2)
}

prob_mat_to_label <- function (xgb, xgb_test, nr, nc) {
  pre_xgb <- matrix(predict(xgb, newdata = xgb_test), nrow = nr,
                    ncol = nc,
                    byrow = TRUE)
  colnames(pre_xgb) <- paste0("G", 1:6)
  pre_xgb_label <- apply(pre_xgb, 1, function(x) {c(0:5)[which(x == max(x))]})
  pre_xgb_label <- factor(pre_xgb_label, levels = 0:5, label = colnames(pre_xgb))
  return(list(pre_xgb, pre_xgb_label))
}

get_pred_summary <- function(res_list, key_names) {
  res_list$summary <- list(rep("", length(res_list$pred[[1]][[2]])))
  for (i in key_names) {
    res_list$summary[[1]] <- paste(res_list$summary[[1]], as.character(res_list$pred[[i]][[2]]), sep = "|")
  }
  res_list$summary <- as.character(lapply(res_list$summary[[1]], function(x) {
    votes <- str_split(x, fixed("|"))[[1]]
    votes <- votes[votes != ""]
    votes <- table(votes)
    paste(names(votes), votes, sep = ":", collapse = "|")
  }))
  return(res_list)
}
