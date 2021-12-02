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
  cluster_genes <- system.file("extdata", "cluster_genes.csv", package = "clinaml")
  cluster_genes <- read.csv(cluster_genes)
  row.names(cluster_genes) <- cluster_genes[,1]

  if (any(is.na(as.numeric(x[,1])))) {
    x <- x[!duplicated(x[,1]),]
    row.names(x) <- x[,1]
    x <- x[,-1]
  }
  if (all(!cluster_genes[,1] %in% colnames(x)) && all(!cluster_genes[,2] %in% colnames(x))) {
    x <- t(x)
  }
  if (any(str_detect(colnames(x), "^ENSG"))) {
    x <- x[,cluster_genes$ens]
    colnames(x) <- cluster_genes[,2]
  } else {
    id_map <- system.file("extdata", "id_map.csv", package = "clinaml")
    id_map <- read.csv(id_map)
    id_map_alias <- str_split(id_map[,4], fixed("|"))
    names(id_map_alias) <- id_map[,3]
    colnames(x) <- str_replace_all(colnames(x), fixed("."), "-")
    idx <- c()
    for (i in 1:ncol(x)) {
      idx <- c(idx, is.character(x[,i]))
    }
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
  }
  x <- x[,cluster_genes$symbol]
  for (i in 1:ncol(x)) {
    x[,i] <- as.numeric(x[,i])
  }
  fil <- apply(x, 2, function(x){all(is.na(x))})
  if (any(fil)) x <- x[,!fil]
  if (max(var(t(x))) > 1000) {
    x <- log2(x + 1)
  }
  return(x)
}
