#! /usr/bin/env Rscript

set.seed(1234)

library(readr)
library(biomaRt)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(RColorBrewer)
library(rstatix)
library(data.table)
library(stringi)
library(stringr)
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(igraph)
library(Matrix)
library(ggplot2)
library(fgsea)
library(BioNet)
library(igraph)
library(msigdbr)
library(dplyr)
library(tibble)
library(future.apply)

dir.create("output")

load(file = "/beegfs/prj/Dewenter_TAC_Backs_lab/Analyses_Gjerga_Questions/all_results.RData") # For this check: https://zenodo.org/records/11208249

########## Prepare databases
## BIOGRID
preview <- readLines("src/BIOGRID-ORGANISM-Mus_musculus-5.0.253.tab.txt", n = 100)
header_line <- grep("^INTERACTOR_A", preview)
biogrid <- fread("src/BIOGRID-ORGANISM-Mus_musculus-5.0.253.tab.txt",
                 skip = header_line - 1)
result <- biogrid[, .(OFFICIAL_SYMBOL_A, OFFICIAL_SYMBOL_B)]
print(head(result))
idx2rem <- which(result$OFFICIAL_SYMBOL_A == result$OFFICIAL_SYMBOL_B)
biogrid <- result[-idx2rem, ]
colnames(biogrid) <- c("Gene1", "Gene2")
biogrid <- unique(biogrid)

biogrid$Gene1 <- str_to_title(biogrid$Gene1)
biogrid$Gene2 <- str_to_title(biogrid$Gene2)

bb <- matrix(data = , nrow = 200000, ncol = 2)
for(ii in 1:nrow(biogrid)){
  bb[ii, ] <- sort(c(biogrid$Gene1[ii], biogrid$Gene2[ii]))
}
bb <- unique(bb)
bb <- as.data.frame(bb)
colnames(bb) <- colnames(biogrid)
biogrid <- bb
biogrid <- biogrid[complete.cases(biogrid), ]


## STRING
string_db <- read.csv("src/10090.protein.physical.links.v12.0.txt", sep="")

ensp <- gsub(pattern = "10090.", replacement = "", x = unique(c(string_db$protein1, string_db$protein2)), fixed = TRUE)

ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl", host="https://www.ensembl.org")
results <- getBM(attributes = c("ensembl_peptide_id", "external_gene_name"),
                 filters = "ensembl_peptide_id",
                 values = ensp,
                 mart = ensembl)

string_db$protein1 <- gsub(pattern = "10090.", replacement = "", x = string_db$protein1, fixed = TRUE)
string_db$protein2 <- gsub(pattern = "10090.", replacement = "", x = string_db$protein2, fixed = TRUE)

mapped_step1 <- merge(string_db, results,
                      by.x = "protein1",
                      by.y = "ensembl_peptide_id",
                      all.x = TRUE)
names(mapped_step1)[names(mapped_step1) == "external_gene_name"] <- "gene1"
final_mapped <- merge(mapped_step1, results,
                      by.x = "protein2",
                      by.y = "ensembl_peptide_id",
                      all.x = TRUE)
names(final_mapped)[names(final_mapped) == "external_gene_name"] <- "gene2"

head(final_mapped)

string_db <- final_mapped[, c(4, 5, 3)]

rm(final_mapped)

string_db <- string_db[, 1:2]


## CAPRINIM
extracted_relations <- read.delim2("src/extracted_relations.txt")
extracted_relations <- unique(extracted_relations)

caprinim <- unique(extracted_relations[, 3:4])

colnames(biogrid) <- c("A", "B")
colnames(string_db) <- c("A", "B")
colnames(caprinim) <- c("A", "B")

caprinim_nodes <- unique(c(caprinim$A, caprinim$B))

idx <- c(intersect(which(string_db$A %in% caprinim_nodes), which(string_db$B %in% caprinim_nodes)),
         intersect(which(string_db$B %in% caprinim_nodes), which(string_db$A %in% caprinim_nodes)))
string_db <- string_db[idx, ]

idx <- c(intersect(which(biogrid$A %in% caprinim_nodes), which(biogrid$B %in% caprinim_nodes)),
         intersect(which(biogrid$B %in% caprinim_nodes), which(biogrid$A %in% caprinim_nodes)))
biogrid <- biogrid[idx, ]

########################################################
## INPUTS
########################################################
ppi_list <- list(
  caprinim = caprinim,
  biogrid = biogrid,
  string  = string_db
)
ppi_list$caprinim <- ppi_list$caprinim[complete.cases(ppi_list$caprinim), ]
ppi_list$biogrid <- ppi_list$biogrid[complete.cases(ppi_list$biogrid), ]
ppi_list$string <- ppi_list$string[complete.cases(ppi_list$string), ]
dim(ppi_list$caprinim)
dim(ppi_list$biogrid)
dim(ppi_list$string)
save(ppi_list, file = "output/ppi_list.RData")
