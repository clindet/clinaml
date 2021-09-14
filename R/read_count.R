#' Function to convert RNA-Seq reads to gene expression matrix (Salmon)
#'
#' @param count_files Character vector for salmon output count files
#' @param samples If count_files is NULL, samples was used to
#' set the value with indir
#' @param indir Input salmon counts dir
#' @param outprefix Output prefix
#' @param txgene_fn Two columns file contains TXNAME (e.g. ENST00000373020.9) and
#' GENEID (ENSG00000000003.15)
#' @param gene2symbol Two columns file contains GENEID (e.g. ENSG00000000003.15) and
#' SYMBOL (TSPAN6)
#' @param workers Number of thread for DESeq2 process
#' @param ... countsFromAbundance ("no", "scaledTPM", "lengthScaledTPM", "dtuScaledTPM")
#' and ignoreAfterBar (TRUE or FALSE). Default is lengthScaledTPM and TRUE
#' @examples
#' \dontrun{
#'   gene2symbol <- fread('~/env/genome/gencode/hg38/v34/tximport.geneid2symbol',
#'     header = T, data.table = FALSE)
#'   txgene_fn <- "~/env/genome/gencode/hg38/v34/tximport.genes"
#'   tx2gene <- fread(txgene_fn)
#'
#'   samples <- readLines("samples")
#'   sampledata <- data.frame(sampleName = samples, condition = rep("disease", length(samples)))
#'   salmon_dir <- "../salmon/output/salmonGTF/"
#'   salmon_files <- file.path(salmon_dir, samples, "quant.sf")
#'   convert_salmon(salmon_files,
#'     outprefix = "salmon/aml",
#'     txgene_fn = txgene_fn,
#'     gene2symbol = gene2symbol
#'   )
#' }
#' @export
convert_salmon <- function (count_files = NULL,
                            samples = NULL,
                            indir = NULL,
                            outprefix = NULL,
                            txgene_fn = "",
                            gene2symbol = "",
                            workers = detectCores(),
                            ...) {
  if (is.null(count_files)) {
    count_files<- file.path(indir, samples, "quant.sf")
  }
  names(count_files) <- samples
  txi <- counts2final(count_files, "salmon", tx2gene, gene2symbol, outprefix, ...)
  txi[[4]] <- deseq2vsd(sampledata, txi[[1]], gene2symbol, paste0(outprefix, "-deseq2vsd.csv"), workers = workers)
  names(txi) <- c("tx.counts", "counts", "tpm", "deseq2vsd")
  return(txi)
}

#' Function to convert RNA-Seq reads to gene expression matrix (Salmon)
#'
#' @param count_files Character vector for salmon output count files
#' @param samples If count_files is NULL, samples was used to
#' set the value with indir
#' @param indir Input salmon counts dir
#' @param outprefix Output prefix
#' @param txgene_fn Two columns file contains TXNAME (e.g. ENST00000373020.9) and
#' GENEID (ENSG00000000003.15)
#' @param gene2symbol Two columns file contains GENEID (e.g. ENSG00000000003.15) and
#' SYMBOL (TSPAN6)
#' @param workers Number of thread for DESeq2 process
#' @param ... countsFromAbundance ("no", "scaledTPM", "lengthScaledTPM", "dtuScaledTPM")
#' and ignoreAfterBar (TRUE or FALSE). Default is lengthScaledTPM and TRUE
#' @examples
#' \dontrun{
#'   gene2symbol <- fread('~/env/genome/gencode/hg38/v34/tximport.geneid2symbol',
#'     header = T, data.table = FALSE)
#'   txgene_fn <- "~/env/genome/gencode/hg38/v34/tximport.genes"
#'   tx2gene <- fread(txgene_fn)
#'
#'   samples <- readLines("samples")
#'   sampledata <- data.frame(sampleName = samples, condition = rep("disease", length(samples)))
#'   salmon_dir <- "../salmon/output/kallistoGTF/"
#'   salmon_files <- file.path(salmon_dir, samples, "abundance.h5")
#'   convert_kallisto(salmon_files,
#'     outprefix = "kallisto/aml",
#'     txgene_fn = txgene_fn,
#'     gene2symbol = gene2symbol
#'   )
#' }
#' @export
convert_kallisto <- function (count_files = NULL,
                              samples = NULL,
                              indir = NULL,
                              outprefix = NULL,
                              txgene_fn = "",
                              gene2symbol = "",
                              workers = detectCores(),
                              ...) {
  if (is.null(count_files)) {
    count_files<- file.path(indir, samples, "abundance.h5")
  }
  names(count_files) <- samples
  txi <- counts2final(count_files, "kallisto", tx2gene, gene2symbol, outprefix, ...)
  txi[[4]] <- deseq2vsd(sampledata, txi[[1]], gene2symbol,
                        paste0(outprefix, "-deseq2vsd.csv"), workers = workers)
  names(txi) <- c("tx.counts", "counts", "tpm", "deseq2vsd")
  return(txi)
}

#' Function to convert RNA-Seq reads to gene expression matrix (featureCounts)
#'
#' @param count_files Character vector for salmon output count files
#' @param samples If count_files is NULL, samples was used to
#' set the value with indir
#' @param indir Input salmon counts dir
#' @param meanFragmentLength Character vector for sequencing library
#' length of samples (e.g. 50, 75, 150, 200)
#' @param out_raw_counts Output path of raw counts table with extra fields
#' @param outprefix Output prefix
#' @param gene2symbol Two columns file contains GENEID (e.g. ENSG00000000003.15) and
#' SYMBOL (TSPAN6)
#' @param workers Number of thread for DESeq2 process
#'
#' @examples
#'
#' samples <- readLines("samples")
#' countfiles <- file.path("../featureCounts", paste0(samples, ".txt"))
#' out_raw_counts <- "featurecounts/raw_counts_584.txt"
#' meanFragmentLength <- rep(150, length(samples))
#' meanFragmentLength[366:477] <- 75
#' gene2symbol <- fread('~/env/genome/gencode/hg38/v34/tximport.geneid2symbol',
#'     header = T, data.table = FALSE)
#' convert_featurecounts(countfiles,
#'   meanFragmentLength = meanFragmentLength,
#'   out_raw_counts = out_raw_counts,
#'   outprefix = "featurecounts/aml-",
#'   gene2symbol = gene2symbol)
convert_featurecounts <- function (countfiles = NULL,
                                   samples = NULL,
                                   indir = NULL,
                                   meanFragmentLength = NULL,
                                   out_raw_counts = "",
                                   outprefix = NULL,
                                   gene2symbol = "",
                                   workers = detectCores()) {
  if (is.null(countfiles)) {
    countfiles <- file.path(indir, paste0(samples, ".txt"))
  }
  if (file.exists(out_raw_counts)) {
    featurecounts <- fread(out_raw_counts, data.table = FALSE)
  } else {
    featurecounts <- featurecounts_preprocess(sampledata, countfiles, out_raw_counts)
  }
  featurecounts <- merge(gene2symbol, featurecounts, by = 1)
  metadata <- featurecounts[,1:7]
  countdata <- featurecounts[,8:ncol(featurecounts)]
  row.names(countdata) <- featurecounts[,1]
  tpm <- as.data.frame(counts2tpm(metadata, countdata, meanFragmentLength))
  filzero <- apply(tpm[,-c(1:2)], 1, function(x) {all(x == 0)})
  fwrite(featurecounts[,-c(3:6)], paste0(outprefix, "-counts.csv"))
  fwrite(tpm[!filzero,], paste0(outprefix, "-tpm.csv"))
  deseq <- deseq2vsd2(sampledata, countdata, gene2symbol,
                      paste0(outprefix, "-deseq2vsd.csv"), workers = workers)
  res <- list(featurecounts[,-c(3:6)], tpm[!filzero,], deseq)
  names(res) <- c("counts", "tpm", "deseq2vsd")
  return(res)
}

#' Function to convert RNA-Seq reads to gene expression matrix (HTseq)
#'
#' @param count_files Character vector for salmon output count files
#' @param samples If count_files is NULL, samples was used to
#' set the value with indir
#' @param indir Input salmon counts dir
#' @param meanFragmentLength Character vector for sequencing library
#' length of samples (e.g. 50, 75, 150, 200)
#' @param exon_length Exon length annotation
#' @param outprefix Output prefix
#' @param gene2symbol Two columns file contains GENEID (e.g. ENSG00000000003.15) and
#' SYMBOL (TSPAN6)
#' @param workers Number of thread for DESeq2 process
#'
#' @examples
#'
#' samples <- readLines("samples")
#' countfiles <- file.path("../featureCounts", paste0(samples, ".txt"))
#' out_raw_counts <- "featurecounts/raw_counts_584.txt"
#' meanFragmentLength <- rep(150, length(samples))
#' meanFragmentLength[366:477] <- 75
#' gene2symbol <- fread('~/env/genome/gencode/hg38/v34/tximport.geneid2symbol',
#'     header = T, data.table = FALSE)
#' convert_htseq(countfiles,
#'   meanFragmentLength = meanFragmentLength,
#'   outprefix = "htseq/aml-",
#'   gene2symbol = gene2symbol)
convert_htseq <- function (countfiles = NULL,
                           samples = NULL,
                           indir = NULL,
                           meanFragmentLength = NULL,
                           exon_length = NULL,
                           outprefix = NULL,
                           gene2symbol = "",
                           workers = detectCores()) {
  htseq <- NULL
  if (is.null(countfiles)) {
    countfiles <- sprintf("../htseq/%s.txt", samples)
  }
  for (i in countfiles) {
    tmp <- fread(i, data.table = FALSE)
    tmp <- tmp[,-1]
    htseq <- cbind(htseq, tmp)
    colnames(htseq)[ncol(htseq)] <- str_remove(i, ".txt$")
  }
  tmp <- fread(i, data.table = FALSE)

  htseq <- cbind(tmp[1], htseq)
  htseq <- merge(exon_length, htseq, by.x = 1, by.y = 1)
  htseq <- merge(gene2symbol, htseq, by.x = 1, by.y = 1)
  metadata <- htseq[,1:3]
  countdata <- htseq[,-c(1:3)]
  row.names(countdata) <- htseq[,1]
  tpm <- as.data.frame(counts2tpm(metadata, countdata, meanFragmentLength))
  filzero <- apply(tpm[,-c(1:2)], 1, function(x) {all(x == 0)})
  fwrite(htseq, paste0(outprefix, "-counts.csv"))
  fwrite(tpm[!filzero,], paste0(outprefix, "-tpm.csv"))
  deseq <- deseq2vsd2(sampledata, countdata, gene2symbol,
                      paste0(outprefix, "-deseq2vsd.csv"))
  res <- list(htseq, tpm[!filzero,], deseq)
  names(res) <- c("counts", "tpm", "deseq2vsd")
  return(res)
}
