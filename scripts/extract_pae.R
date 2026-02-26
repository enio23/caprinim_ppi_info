library(readr)
library(jsonlite)
library(purrr)

dd <- list.dirs(path = "resource_chunk_1/")
dd <- dd[which(grepl(pattern = "pos_", x = dd, fixed = TRUE))]

dir.create("PAE_Evals")


#### Find PAE
# ---- collect ALL values for any of these keys anywhere in the JSON ----
find_all_keys <- function(x, targets) {
  out <- list()
  walk <- function(node) {
    if (is.list(node)) {
      if (!is.null(names(node))) {
        hit <- intersect(names(node), targets)
        if (length(hit) > 0) {
          for (k in hit) out[[length(out) + 1]] <<- list(key = k, value = node[[k]])
        }
      }
      for (v in node) walk(v)
    }
  }
  walk(x)
  out
}

# ---- does this look like a square numeric "list-of-rows"? ----
looks_like_square_list <- function(v) {
  if (!is.list(v) || length(v) < 2) return(FALSE)
  lens <- vapply(v, length, integer(1))
  if (length(unique(lens)) != 1) return(FALSE)
  n <- length(v)
  if (lens[1] != n) return(FALSE)
  # check that entries can be numeric
  ok <- vapply(v, function(row) {
    suppressWarnings(all(is.finite(as.numeric(row)) | is.na(as.numeric(row))))
  }, logical(1))
  all(ok)
}

# ---- coerce candidate into numeric matrix if possible ----
coerce_to_numeric_matrix <- function(v) {
  # already matrix
  if (is.matrix(v)) {
    m <- v
    storage.mode(m) <- "double"
    return(m)
  }
  # data.frame
  if (is.data.frame(v)) {
    m <- as.matrix(v)
    storage.mode(m) <- "double"
    return(m)
  }
  # list-of-rows
  if (looks_like_square_list(v)) {
    m <- do.call(rbind, lapply(v, function(r) as.numeric(r)))
    storage.mode(m) <- "double"
    return(m)
  }
  # sometimes it's a 1-element wrapper list
  if (is.list(v) && length(v) == 1) {
    return(coerce_to_numeric_matrix(v[[1]]))
  }
  stop("Candidate is not a numeric NxN matrix-like structure.")
}

# ---- v1 flattened format {residue1,residue2,distance} ----
pae_from_v1 <- function(obj) {
  r1 <- as.integer(obj$residue1)
  r2 <- as.integer(obj$residue2)
  d  <- as.numeric(obj$distance)
  n  <- max(c(r1, r2), na.rm = TRUE)
  m  <- matrix(NA_real_, nrow = n, ncol = n)
  m[cbind(r1, r2)] <- d
  m
}

extract_pae_from_full_data <- function(json_path) {
  jd <- fromJSON(json_path, simplifyVector = FALSE)
  
  # ---- 1) try v2-like matrix keys, but pick the first *matrix-like* candidate ----
  targets_v2 <- c("pae", "predicted_alignment_error", "predicted_aligned_error")
  cand <- find_all_keys(jd, targets_v2)
  
  if (length(cand) > 0) {
    for (i in seq_along(cand)) {
      m <- tryCatch(coerce_to_numeric_matrix(cand[[i]]$value), error = function(e) NULL)
      if (!is.null(m) && nrow(m) == ncol(m) && nrow(m) > 1) {
        return(list(pae = m, format = "v2", key = cand[[i]]$key))
      }
    }
  }
  
  # ---- 2) try v1 flattened keys anywhere ----
  r1 <- find_all_keys(jd, c("residue1"))
  r2 <- find_all_keys(jd, c("residue2"))
  dd <- find_all_keys(jd, c("distance"))
  
  if (length(r1) > 0 && length(r2) > 0 && length(dd) > 0) {
    # use the first triple that gives a square-ish matrix
    obj <- list(residue1 = r1[[1]]$value, residue2 = r2[[1]]$value, distance = dd[[1]]$value)
    m <- pae_from_v1(obj)
    return(list(pae = m, format = "v1", key = "residue1/residue2/distance"))
  }
  
  stop("No usable PAE found in ", json_path,
       "\nTip: run str(fromJSON(json_path, simplifyVector=FALSE), max.level=3) to inspect.")
}

for(ii in 1:length(dd)){
  
  print(paste0("Step ---- ", ii, "/", length(dd)))
  
  # base_dir <- gsub(pattern = "src//", replacement = "", x = dd[ii], fixed = TRUE)
  
  # dir.create(paste0("PAE_Evals/", base_dir))
  
  ff <- list.files(dd[ii])
  ff <- ff[which(grepl(pattern = "full_data", x = ff, fixed = TRUE))]
  base_file <- gsub(pattern = ".json", replacement = "", x = ff, fixed = TRUE)
  
  res <- extract_pae_from_full_data(paste0(dd[ii], "/", ff))
  
  pae <- res$pae
  
  write.csv(pae, paste0("PAE_Evals/", base_file, "_pae.csv"), row.names = FALSE)
  
}