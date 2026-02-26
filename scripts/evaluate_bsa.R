#! /usr/bin/env Rscript

library(bio3d)
library(readr)

dir.create("BSA_Evals")

aa3_to1 <- function(x, unknown = "X") {
  # mapping (case-insensitive keys)
  map <- c(
    Ala="A", Arg="R", Asn="N", Asp="D", Cys="C",
    Gln="Q", Glu="E", Gly="G", His="H", Ile="I",
    Leu="L", Lys="K", Met="M", Phe="F", Pro="P",
    Ser="S", Thr="T", Trp="W", Tyr="Y", Val="V",
    Sec="U", Pyl="O", Asx="B", Glx="Z", Xle="J",
    Ter="*", Stop="*", Xxx="X"
  )
  
  # common PDB variants → canonical 3-letter
  # (add more as needed)
  synonyms <- c(
    HID="His", HIE="His", HIP="His", HSD="His", HSE="His",
    MSE="Met", SEL="Sec", TPO="Thr", SEP="Ser", PTR="Tyr",
    HYP="Pro"
  )
  
  # Normalize one input (string or vector of tokens) -> vector of 3-letter tokens
  normalize_tokens <- function(v) {
    if (length(v) == 1L) {
      # split a chain like "Met-Gly Asp" -> c("Met","Gly","Asp")
      toks <- unlist(strsplit(v, "[^A-Za-z]+"))
      toks <- toks[nzchar(toks)]
    } else {
      toks <- v
    }
    toks <- sub("^([A-Za-z])(.*)$", "\\U\\1\\L\\2", toks, perl = TRUE) # Title Case
    toks
  }
  
  to1 <- function(one) {
    toks <- normalize_tokens(one)
    
    # apply synonyms first
    canon <- ifelse(tolower(toks) %in% tolower(names(synonyms)),
                    synonyms[match(tolower(toks), tolower(names(synonyms)))],
                    toks)
    
    # map to one-letter
    ol <- ifelse(tolower(canon) %in% tolower(names(map)),
                 map[match(tolower(canon), tolower(names(map)))],
                 unknown)
    
    paste0(ol, collapse = "")
  }
  
  # vectorized: each element of x becomes one output sequence
  vapply(x, to1, FUN.VALUE = character(1))
}

dirs <- list.dirs(path = "resource_chunk_1/")
dirs <- dirs[2:length(dirs)]
dirs <- gsub(pattern = "//", replacement = "/", x = dirs) 

mouse_ids <- gsub(pattern = "resource_chunk_1//pos_", replacement = "", x = dirs, fixed = TRUE)
mouse_ids <- toupper(unique(unlist(strsplit(x = mouse_ids, split = "_", fixed = TRUE))))

mart <- biomaRt::useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")

attrs <- c(
  "ensembl_transcript_id","uniprotswissprot","uniprotsptrembl",
  "pfam","pfam_start","pfam_end",
  "ensembl_gene_id","external_gene_name","external_synonym","description"
)

ensembl_df <- biomaRt::getBM(attributes = attrs, filters = "uniprotswissprot", values = mouse_ids, mart = mart)


source("compute_bsa.R")

for(ii in 1:length(dirs)){
  
  print(paste0("Step ---- ", ii, "/", length(dirs)))
  
  case <- strsplit(x = dirs[ii], split = "/", fixed = TRUE)[[1]][4]
  
  A = toupper(strsplit(x = case, split = "_", fixed = TRUE)[[1]][2])
  B = toupper(strsplit(x = case, split = "_", fixed = TRUE)[[1]][3])
  
  domainsA <- unique(ensembl_df[which(ensembl_df$uniprotswissprot == A), 4:6])
  domainsB <- unique(ensembl_df[which(ensembl_df$uniprotswissprot == B), 4:6])
  domain_map <- NULL
  if(nrow(domainsA) == 0){
    if(nrow(domainsB) > 0){
      domainsB$chain <- "B"
      domain_map <- domainsB
    }
  }
  if(nrow(domainsB) == 0){
    if(nrow(domainsA) > 0){
      domainsA$chain <- "A"
      domain_map <- domainsA
    }
  }
  if((nrow(domainsA)) > 0 && (nrow(domainsB) > 0)){
    domainsA$chain <- "A"
    domainsB$chain <- "B"
    domain_map <- rbind(domainsA, domainsB)
  }
  
  pdb_file <- list.files(path = dirs[ii])
  pdb_path <- paste0(dirs[ii], "/", pdb_file[which(grepl(pattern = ".pdb", x = pdb_file, fixed = TRUE))[1]])
  
  if(file.exists(pdb_path)){
    
    res <- compute_bsa(pdb_path = pdb_path, domain_map = domain_map, cleanup = TRUE)
    
    res$aa <- aa3_to1(x = res$chain)
    
    rel <- gsub(pattern = "resource_chunk_1/pos_", replacement = "", x = dirs[ii], fixed = TRUE)
    
    write_csv(x = res, file = paste0("BSA_Evals/", rel, ".csv"))
    
  }
  
}

