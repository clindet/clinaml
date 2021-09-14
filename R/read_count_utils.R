featurecounts_preprocess <- function (sampledata, countfiles, outfn) {
  final <- NULL
  for (i in 1:length(countfiles)) {
    print(i)
    tmp <- as.data.frame(data.table::fread(countfiles[i]))
    colnames(tmp)[7] <- sampledata$sampleName[i]
    if (is.null(final)) {
      final <- tmp
    } else {
      final <- merge(final, tmp, by=1:6)
    }
  }
  fwrite(final, outfn, sep = '\t', row.names = F, quote = F)
  return(final)
}

counts2tpm <- function(metadata, counts, meanFragmentLength) {
  # Ensure valid arguments.
  stopifnot(length(metadata$Length) == nrow(counts))
  stopifnot(length(meanFragmentLength) == ncol(counts))

  # Compute effective lengths of features in each library.
  effLen <- do.call(cbind, lapply(1:ncol(counts), function(i) {
    metadata$Length - meanFragmentLength[i] + 1
  }))

  # Exclude genes with length less than the mean fragment length.
  idx <- apply(effLen, 1, function(x) min(x) > 1)
  counts <- counts[idx,]
  effLen <- effLen[idx,]
  metadata <- metadata[idx,]

  # Process one column at a time.
  tpm <- do.call(cbind, lapply(1:ncol(counts), function(i) {
    rate = log(counts[,i]) - log(effLen[,i])
    denom = log(sum(exp(rate)))
    exp(rate - denom + log(1e6))
  }))

  # Copy the column names from the original matrix.
  colnames(tpm) <- colnames(counts)
  tpm <- data.frame(metadata[,1:2], tpm)
  return(tpm)
}


deseq2vsd <- function (sampledata, txi, gene2symbol, outfn, workers = 39) {
  ddsTxi <- DESeqDataSetFromTximport(txi, sampledata, ~1)
  BPPARAM <- MulticoreParam(workers = workers)
  dds <- DESeq(ddsTxi, parallel = TRUE, BPPARAM = BPPARAM)
  vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
  noZero <- (rowSums(counts(dds)) > 0)
  dat.out <- assay(vsd[noZero, ])
  dat.out <- merge(gene2symbol, dat.out, by.x = 1, by.y = 0, all.y = T)
  saveRDS(txi, paste0(outfn, ".txi.rds"))
  fwrite(dat.out, file = outfn)
  return(dat.out)
}

deseq2vsd2 <- function (sampledata, countdata, gene2symbol, outfn, workers = 39) {
  ddsHTSeq <- DESeqDataSetFromMatrix(countdata, colData = sampledata, design=~1)
  BPPARAM <- MulticoreParam(workers = workers)
  dds <- DESeq(ddsHTSeq, parallel = TRUE, BPPARAM = BPPARAM)
  fil <- rowSums(counts(dds)) > 0
  dds <- dds[fil,]
  vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
  noZero <- (rowSums(counts(dds))>0)
  dat.out <- assay(vsd[noZero,])
  dat.out <- merge(gene2symbol, dat.out, by.x = 1, by.y = 0, all.y = T)
  fwrite(dat.out, file = outfn)
  return(dat.out)
}

counts2final <- function (count.files, type, tx2gene, gene2symbol, outprefix,
                          countsFromAbundance = "lengthScaledTPM",
                          ignoreAfterBar = TRUE) {
  tx.counts <- tximport(count.files, type = type, tx2gene = tx2gene,
                        countsFromAbundance = countsFromAbundance,
                        ignoreAfterBar = ignoreAfterBar)
  counts <- merge(gene2symbol, tx.counts$counts, by.x = 1, by.y = 0)
  tpm <- merge(gene2symbol, tx.counts$abundance, by.x = 1, by.y = 0)
  filzero <- apply(tpm[,-c(1:2)], 1, function(x) {all(x == 0)})
  if (!is.null(outprefix)) {
    fwrite(counts, paste0(outprefix, "-counts.csv"))
    fwrite(tpm[!filzero,], paste0(outprefix, "-tpm.csv"))
  }
  return(list(tx.counts, counts, tpm[!filzero,]))
}
