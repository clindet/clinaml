#' Function to convert GTF to tximport.geneid2symbol and tximport.genes
#'
#' @param gtf GTF path (e.g. GENECODE v34)
#' @param export
#' @examples
#' \dontrun {
#'   gtf2tximport("gencode.v38.annotation.gtf")
#'   gtf2tximport("tximport.genes", use_remote_genecode = TRUE)
#' }
gtf2tximport <- function (gtf = NULL, outprefix = "tximport",
                          use_remote_genecode = FALSE, version = 38) {
  gtf <- read_gtf(gtf, use_remote_genecode = use_remote_genecode, version = version)
  symbol <- str_split(gtf[gtf[,3] == "gene", 9], "; ")
  GENEID <- sapply(symbol, function(x) {return(x[1])})
  SYMBOL <- sapply(symbol, function(x) {return(x[3])})
  GENEID <- str_remove_all(GENEID, '"| |gene_id')
  SYMBOL <- str_remove_all(SYMBOL, '"| |gene_name')
  fwrite(data.frame(GENEID, SYMBOL), paste0(outprefix, ".geneid2symbol"),
         sep = ",", quote = FALSE, row.names = FALSE)

  gtf <- str_split(gtf[gtf[,3] %in% c("transcript",
                                      "mRNA"), 9], "; ")
  GENEID <- sapply(gtf, function(x) {x[1]})
  TXNAME <- sapply(gtf, function(x) {x[2]})
  GENEID <- str_remove_all(GENEID, '"| |gene_id')
  TXNAME <- str_remove_all(TXNAME, '"| |transcript_id')
  fwrite(data.frame(TXNAME, GENEID), paste0(outprefix, ".genes"),
         sep = ",", quote = FALSE, row.names = FALSE)
}

#' Function to convert GTF to exon length
#'
#' @param gtf GTF path (e.g. GENECODE v34)
gtf2exon_len <- function (gtf = NULL, outprefix = "exon",
                          use_remote_genecode = FALSE, version = 38) {
  gtf <- read_gtf(gtf, use_remote_genecode = use_remote_genecode, version = version)
}

read_gtf <- function (gtf = NULL,
                      use_remote_genecode = FALSE, version = 38) {
  url <- NULL
  if (str_detect(gtf, "http:|https:|ftp:")) {
    url <- gtf
  }
  if (use_remote_genecode) {
    url <- sprintf("http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_%s/gencode.v%s.annotation.gtf.gz",
                   version, version)
  }
  if (!is.null(url)) {
    download.file(url, basename(url))
    gtf <- basename(url)
  }
  if (str_detect(gtf, ".gz$")) {
    gunzip(gtf)
    gtf <- str_remove(gtf, ".gz$")
  }
  gtf <- suppressWarnings(fread(gtf, data.table = FALSE))
}
