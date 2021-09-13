library(stringr)
library(clusterProfiler)
library(org.Hs.eg.db)

exp <- read.csv("tcga_deseq2_test.csv", check.names = FALSE)
dat <- read.csv("cluster_genes.csv")

dat <- dat[dat[,2] %in% colnames(exp),]


dat$ens <- str_remove_all(dat$ens, "[.][0-9]*")

map_dt <- bitr(dat$ens, fromType = "ENSEMBL",toType = c("ENTREZID"), OrgDb = org.Hs.eg.db)

dat <- merge(dat, map_dt, by = 1)

map_dt <- bitr(dat$ens, fromType = "ENSEMBL",toType = c("ALIAS"), OrgDb = org.Hs.eg.db)

dat$alias <- NA
for (i in 1:nrow(dat)) {
  idx <- map_dt[,1] == dat$ens[i]
  dat$alias[i] <- paste(map_dt[idx,2], collapse = "|")
}

write.csv(dat, "id_map.csv", row.names = FALSE, quote = FALSE)
