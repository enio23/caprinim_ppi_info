#! /usr/bin/env Rscript

library(bio3d)
library(readr)
library(jsonlite)

dir.create("AF_Features")

dirs <- list.dirs(path = "resource_chunk_12/")
dirs <- dirs[2:length(dirs)]
dirs <- gsub(pattern = "//", replacement = "/", x = dirs) 

af_features <- matrix(data = , nrow = length(dirs), ncol = 5)
for(ii in 1:length(dirs)){
  
  print(paste0("Step ---- ", ii, "/", length(dirs)))
  
  ff <- list.files(path = paste0(dirs[ii], "/"))
  ff <- ff[which(grepl(pattern = "summary_confidences", x = ff, fixed = TRUE))]
  summary_confidences <- fromJSON(paste0(dirs[ii], "/", ff))
  
  # Features
  chain_ptm_a <- summary_confidences$chain_ptm[1]
  chain_ptm_b <- summary_confidences$chain_ptm[2]
  min_pae_ab <- summary_confidences$chain_pair_pae_min[1, 2]
  ptm <- summary_confidences$ptm
  
  ff <- gsub(pattern = "_summary_confidences_0.json", replacement = "", x = ff, fixed = TRUE)
  
  af_features[ii, ] <- c(gsub(pattern = "fold_", replacement = "", x = ff, fixed = TRUE),
                         chain_ptm_a,
                         chain_ptm_b,
                         min_pae_ab,
                         ptm)
  
}

colnames(af_features) <- c("pair_id", "chain_ptm_a", "chain_ptm_b", "min_pae_ab", "ptm")

write.table(x = af_features, file = paste0("AF_Features/resource_chunk_12.txt"), 
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)



#### Merge
tt1 <- read.delim(file = "AF_Features/resource_chunk_1.txt")
tt2 <- read.delim(file = "AF_Features/resource_chunk_2.txt")
tt3 <- read.delim(file = "AF_Features/resource_chunk_3.txt")
tt4 <- read.delim(file = "AF_Features/resource_chunk_4.txt")
tt5 <- read.delim(file = "AF_Features/resource_chunk_5.txt")
tt6 <- read.delim(file = "AF_Features/resource_chunk_6.txt")
tt7 <- read.delim(file = "AF_Features/resource_chunk_7.txt")
tt8 <- read.delim(file = "AF_Features/resource_chunk_8.txt")
tt9 <- read.delim(file = "AF_Features/resource_chunk_9.txt")
tt10 <- read.delim(file = "AF_Features/resource_chunk_10.txt")
tt11 <- read.delim(file = "AF_Features/resource_chunk_11.txt")
tt12 <- read.delim(file = "AF_Features/resource_chunk_12.txt")

tt <- unique(rbind(tt1, tt2, tt3, tt4, tt5, tt6, tt7, tt8, tt9, tt10, tt11, tt12))

write.table(x = tt, file = "AF_Features/AF_Server_Features.txt", quote = FALSE, 
            sep = "\t", row.names = FALSE, col.names = TRUE)

