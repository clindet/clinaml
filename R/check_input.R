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
  id_map <- system.file("extdata", "id_map.csv", package = "clinaml")
  id_map <- read.csv(id_map)
  id_map_alias <- str_split(id_map[,4], fixed("|"))
  names(id_map_alias) <- id_map[,3]
  colnames(x) <- str_replace_all(colnames(x), fixed("."), "-")
  idx <- c()
  for (i in 1:ncol(x)) {
    idx <- c(idx, is.character(x[,i]))
  }
  x <- x[,!idx]
  fil <- str_detect(colnames(x), "^ENSG")
  if (any(fil)) {
    colnames(x)[fil] <- id_map[match(colnames(x)[fil], id_map[,2]),3]
  }
  if (all(str_detect(colnames(x), "^X"))) {
    idx <- match(colnames(x), paste0("X", id_map[,1]))
    colnames(x) <- id_map[idx,3]
  }
  if (!all(colnames(x) %in% id_map[,3])) {
    idx <- which(!colnames(x) %in% id_map[,3])
    for (i in idx) {
      res <- sapply(id_map_alias, function(x2) {
        return(colnames(x)[i] %in% x2)
      })
      if (any(res)) {
        colnames(x)[i] <- unique(names(id_map_alias)[res])
      }
    }
  }
  return(x)
}
