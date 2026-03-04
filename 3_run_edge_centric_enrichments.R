library(tidyverse)
library(glue)
library(stringr)
library(purrr)
library(forcats)
library(ggplot2)
library(graphite)
library(graph)
library(patchwork)

bionet_thresholds <- c(0.05)
graphite_species <- "mmusculus"

collapse_pathway_direction <- TRUE
restrict_to_cardio_pathways <- TRUE

min_pathway_nodes_in_benchmark <- 4
min_pathway_edges_in_universe <- 3

alpha <- 0.1

run_empirical_permutations <- FALSE
n_permutations <- 1000L
set.seed(1)

outdir <- "edge_topology_enrichment_graphite_all_topology_dbs_threshold_0p05"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "plots"), recursive = TRUE, showWarnings = FALSE)
dir.create("cache", showWarnings = FALSE)

# Formatting / plotting safety
fmt_thr <- function(x) formatC(x, format = "f", digits = 2)

empty_plot <- function(title = "No data available", subtitle = NULL) {
  ggplot() +
    annotate("text", x = 0, y = 0, label = "No data available for this panel", size = 5, color = "grey25") +
    xlim(-1, 1) + ylim(-1, 1) +
    labs(title = title, subtitle = subtitle) +
    theme_void(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, color = "grey35")
    )
}

safe_ggsave <- function(filename, plot, width, height, dpi = 320, bg = "white") {
  tryCatch({
    ggsave(filename, plot = plot, width = width, height = height, dpi = dpi, bg = bg)
  }, error = function(e) {
    message(glue("ggsave failed for {filename}: {conditionMessage(e)}"))
    fallback <- empty_plot(
      title = basename(filename),
      subtitle = paste("Plot rendering failed:", conditionMessage(e))
    )
    ggsave(filename, plot = fallback, width = width, height = height, dpi = dpi, bg = bg)
  })
}

sig_stars <- function(p) {
  case_when(
    is.na(p) ~ "ns",
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE ~ "ns"
  )
}

cluster_heatmap_order <- function(df, row_col = "pathway_label", col_col = "resource", value_col = "neglog10_padj") {
  if (nrow(df) == 0) return(list(row_order = character(0), col_order = character(0)))
  
  mat_df <- df %>%
    select(all_of(c(row_col, col_col, value_col))) %>%
    distinct() %>%
    pivot_wider(names_from = all_of(col_col), values_from = all_of(value_col), values_fill = 0)
  
  if (nrow(mat_df) == 0) return(list(row_order = character(0), col_order = character(0)))
  
  rn <- mat_df[[row_col]]
  mat <- mat_df %>% select(-all_of(row_col)) %>% as.data.frame()
  rownames(mat) <- rn
  mat <- as.matrix(mat)
  mat[!is.finite(mat)] <- 0
  
  row_order <- rownames(mat)
  col_order <- colnames(mat)
  
  if (nrow(mat) >= 2) {
    row_hc <- hclust(dist(mat))
    row_order <- rownames(mat)[row_hc$order]
  }
  if (ncol(mat) >= 2) {
    col_hc <- hclust(dist(t(mat)))
    col_order <- colnames(mat)[col_hc$order]
  }
  
  list(row_order = row_order, col_order = col_order)
}

# LOAD REQUIRED OBJECTS
load(file = "output/ppi_list.RData")

if (!exists("ppi_list")) stop("Missing required object: ppi_list")
if (!all(c("caprinim", "biogrid", "string") %in% names(ppi_list))) {
  stop("ppi_list must contain: ppi_list$caprinim, ppi_list$biogrid, ppi_list$string")
}

# Edge canonicalization
canonicalize_edges_df <- function(df, col1, col2, resource_name = NA_character_) {
  stopifnot(col1 %in% names(df), col2 %in% names(df))
  
  df %>%
    transmute(
      gene1_raw = as.character(.data[[col1]]),
      gene2_raw = as.character(.data[[col2]])
    ) %>%
    mutate(
      gene1_raw = trimws(gene1_raw),
      gene2_raw = trimws(gene2_raw)
    ) %>%
    filter(!is.na(gene1_raw), !is.na(gene2_raw), gene1_raw != "", gene2_raw != "") %>%
    filter(gene1_raw != gene2_raw) %>%
    mutate(
      g1u = toupper(gene1_raw),
      g2u = toupper(gene2_raw),
      a = pmin(g1u, g2u),
      b = pmax(g1u, g2u),
      edge_id = paste(a, b, sep = "||")
    ) %>%
    distinct(edge_id, .keep_all = TRUE) %>%
    transmute(gene1 = a, gene2 = b, edge_id = edge_id, resource = resource_name)
}

make_edge_id <- function(a, b, undirected = TRUE) {
  a <- toupper(trimws(as.character(a)))
  b <- toupper(trimws(as.character(b)))
  if (undirected) {
    x <- pmin(a, b)
    y <- pmax(a, b)
    paste(x, y, sep = "||")
  } else {
    paste(a, b, sep = "->")
  }
}

pretty_pathway_name <- function(x) {
  x %>% str_replace_all("_", " ") %>% str_squish()
}

# Cardio pathway filtering (title matching; unchanged)
is_cardio_relevant_pathway <- function(pathway_name) {
  nm <- toupper(pathway_name)
  
  patterns <- c(
    "MUSCLE CONTRACTION",
    "CARDIAC CONTRACTION",
    "CARDIOMYOPATHY",
    "HYPERTROPH",
    "CARDIOMYOCYTE",
    "CONTRACTIL",
    "SARCOMER",
    "MYOFILAMENT",
    "EXCITATION[- ]CONTRACTION",
    "CALCIUM",
    "ARRHYTHM",
    "VENTRICULAR",
    "DILAT"
  )
  
  any(str_detect(nm, regex(paste(patterns, collapse = "|"), ignore_case = TRUE)))
}

# Discover & load as many graphite pathway DBs as possible
`%||%` <- function(a, b) if (!is.null(a)) a else b

discover_graphite_databases <- function(verbose = TRUE) {
  # Try official helper(s) if present; otherwise fallback list.
  fallback <- c(
    "kegg", "reactome", "wikipathways", "panther", "smpdb", "pharmgkb", "pathbank",
    # Older / sometimes-available names (safe: will be tried and skipped on error)
    "biocarta", "nci", "pid", "humancyc", "spike", "netpath", "signalink"
  ) %>% unique()
  
  dbs <- NULL
  
  # Try a few likely function names across graphite versions
  for (fn in c("pathwayDatabases", "pathwayDatabase", "pathwayDBs", "databases")) {
    if (exists(fn, where = asNamespace("graphite"), inherits = FALSE)) {
      obj <- tryCatch(get(fn, envir = asNamespace("graphite"))(), error = function(e) NULL)
      if (!is.null(obj)) {
        if (is.character(obj)) {
          dbs <- obj
          break
        } else if (is.data.frame(obj)) {
          cand <- NULL
          for (col in c("database", "db", "name", "source")) {
            if (col %in% names(obj)) { cand <- obj[[col]]; break }
          }
          dbs <- cand %||% obj[[1]]
          break
        } else if (is.list(obj)) {
          # sometimes list of db names / records
          if (is.character(obj)) {
            dbs <- obj
            break
          }
        }
      }
    }
  }
  
  dbs <- dbs %||% fallback
  dbs <- dbs %>% as.character() %>% tolower() %>% unique() %>% sort()
  
  if (verbose) message(glue("Candidate graphite pathway DBs to try ({length(dbs)}): {paste(dbs, collapse=', ')}"))
  dbs
}

load_graphite_topology_edges_db <- function(species = "mmusculus",
                                            pathway_db = "kegg",
                                            collapse_direction = TRUE,
                                            restrict_cardio = TRUE,
                                            cache_file = NULL,
                                            use_cache = TRUE,
                                            ncpus = max(1, parallel::detectCores() - 1),
                                            verbose = TRUE) {
  pathway_db <- tolower(pathway_db)
  
  if (!is.null(cache_file) && use_cache && file.exists(cache_file)) {
    if (verbose) message(glue("Loading cached topology edges from: {cache_file}"))
    return(readRDS(cache_file))
  }
  
  old_ncpus <- getOption("Ncpus")
  options(Ncpus = ncpus)
  on.exit(options(Ncpus = old_ncpus), add = TRUE)
  
  if (verbose) message(glue("Loading {toupper(pathway_db)} pathways from graphite ({species})..."))
  pw <- graphite::pathways(species, pathway_db)
  
  pw_sym <- tryCatch(
    graphite::convertIdentifiers(pw, "symbol"),
    error = function(e) {
      warning(glue("convertIdentifiers(..., 'symbol') failed for {toupper(pathway_db)}: {conditionMessage(e)}"))
      pw
    }
  )
  
  if (isTRUE(restrict_cardio)) {
    pw_names <- vapply(pw_sym, function(x) x@title, FUN.VALUE = character(1))
    keep <- vapply(pw_names, is_cardio_relevant_pathway, FUN.VALUE = logical(1))
    pw_sym <- pw_sym[keep]
    if (verbose) message(glue("  Retained {sum(keep)} cardio-relevant {toupper(pathway_db)} pathways"))
  } else {
    if (verbose) message(glue("  Using all {length(pw_sym)} {toupper(pathway_db)} pathways"))
  }
  
  pw_all <- as.list(pw_sym)
  if (length(pw_all) == 0) {
    return(tibble(
      source = character(0), target = character(0), edge_id = character(0),
      pathway_id = character(0), pathway_name = character(0), source_db = character(0)
    ))
  }
  
  if (verbose) message(glue("Extracting topology edges from {length(pw_all)} {toupper(pathway_db)} pathways..."))
  
  extract_one_pathway_edges <- function(pw) {
    pathway_id <- tryCatch(as.character(pw@id), error = function(e) NA_character_)
    pathway_name <- tryCatch(as.character(pw@title), error = function(e) NA_character_)
    
    edf <- tryCatch(graphite::edges(pw, which = "proteins"), error = function(e) NULL)
    if (is.null(edf) || nrow(edf) == 0) return(NULL)
    
    nms_low <- tolower(names(edf))
    from_idx <- which(nms_low %in% c("src", "source", "from"))
    to_idx   <- which(nms_low %in% c("dest", "target", "to"))
    
    if (length(from_idx) == 0 || length(to_idx) == 0) {
      if (ncol(edf) < 2) return(NULL)
      from_idx <- 1
      to_idx <- 2
    } else {
      from_idx <- from_idx[1]
      to_idx <- to_idx[1]
    }
    
    tibble(
      source = as.character(edf[[from_idx]]),
      target = as.character(edf[[to_idx]])
    ) %>%
      mutate(
        source = toupper(trimws(source)),
        target = toupper(trimws(target))
      ) %>%
      filter(!is.na(source), !is.na(target), source != "", target != "", source != target) %>%
      mutate(
        edge_id = make_edge_id(source, target, undirected = collapse_direction),
        pathway_id = pathway_id,
        pathway_name = pathway_name,
        source_db = pathway_db
      ) %>%
      distinct(pathway_id, pathway_name, source_db, edge_id, .keep_all = TRUE)
  }
  
  edge_list <- vector("list", length(pw_all))
  for (i in seq_along(pw_all)) {
    if (verbose && (i %% 50 == 0 || i == 1 || i == length(pw_all))) {
      message(glue("  ... pathway {i}/{length(pw_all)}"))
    }
    edge_list[[i]] <- extract_one_pathway_edges(pw_all[[i]])
  }
  
  edge_tbl <- bind_rows(edge_list)
  if (nrow(edge_tbl) == 0) {
    edge_tbl <- tibble(
      source = character(0), target = character(0), edge_id = character(0),
      pathway_id = character(0), pathway_name = character(0), source_db = character(0)
    )
  }
  
  if (!is.null(cache_file)) {
    dir.create(dirname(cache_file), showWarnings = FALSE, recursive = TRUE)
    saveRDS(edge_tbl, cache_file)
    if (verbose) message(glue("Saved cache: {cache_file}"))
  }
  
  edge_tbl
}

build_pathway_topology_reference <- function(pathway_edges_tbl, benchmark_genes, edge_universe_ids,
                                             min_nodes = 4, min_edges = 3) {
  benchmark_genes <- unique(toupper(benchmark_genes))
  edge_universe_ids <- unique(edge_universe_ids)
  
  pathway_ref <- pathway_edges_tbl %>%
    mutate(source = toupper(source), target = toupper(target)) %>%
    filter(source %in% benchmark_genes, target %in% benchmark_genes) %>%
    group_by(source_db, pathway_id, pathway_name) %>%
    summarise(
      pathway_nodes = list(sort(unique(c(source, target)))),
      pathway_edge_ids_topology = list(unique(edge_id)),
      .groups = "drop"
    ) %>%
    mutate(
      n_nodes_in_benchmark = map_int(pathway_nodes, length),
      pathway_edge_ids_in_universe = map(pathway_edge_ids_topology, ~ intersect(.x, edge_universe_ids)),
      n_pathway_edges_in_universe = map_int(pathway_edge_ids_in_universe, length)
    ) %>%
    filter(
      n_nodes_in_benchmark >= min_nodes,
      n_pathway_edges_in_universe >= min_edges
    )
  
  message(glue("Pathways passing benchmark overlap filters: {nrow(pathway_ref)}"))
  pathway_ref
}

# Enrichment statistics
fisher_edge_enrichment <- function(module_edge_ids, edge_universe_ids, pathway_edge_ids) {
  module_edge_ids  <- intersect(unique(module_edge_ids), unique(edge_universe_ids))
  pathway_edge_ids <- intersect(unique(pathway_edge_ids), unique(edge_universe_ids))
  edge_universe_ids <- unique(edge_universe_ids)
  
  U <- length(edge_universe_ids)
  M <- length(module_edge_ids)
  P <- length(pathway_edge_ids)
  X <- length(intersect(module_edge_ids, pathway_edge_ids))
  
  a <- X
  b <- M - X
  c <- P - X
  d <- U - M - P + X
  
  if (any(c(a, b, c, d) < 0) || any(is.na(c(a, b, c, d)))) {
    return(tibble(
      overlap_edges = X, module_edges = M, pathway_edges = P, universe_edges = U,
      odds_ratio = NA_real_, p_value = NA_real_
    ))
  }
  
  ft <- suppressWarnings(
    fisher.test(matrix(c(a, b, c, d), nrow = 2, byrow = TRUE), alternative = "greater")
  )
  
  tibble(
    overlap_edges = X,
    module_edges = M,
    pathway_edges = P,
    universe_edges = U,
    odds_ratio = unname(ft$estimate),
    p_value = ft$p.value
  )
}

empirical_edge_p <- function(module_edge_ids, edge_universe_ids, pathway_edge_ids, n_perm = 1000L) {
  module_edge_ids  <- intersect(unique(module_edge_ids), unique(edge_universe_ids))
  pathway_edge_ids <- intersect(unique(pathway_edge_ids), unique(edge_universe_ids))
  edge_universe_ids <- unique(edge_universe_ids)
  
  m <- length(module_edge_ids)
  if (m == 0L || length(pathway_edge_ids) == 0L || m > length(edge_universe_ids)) return(NA_real_)
  
  obs <- length(intersect(module_edge_ids, pathway_edge_ids))
  
  hits <- replicate(n_perm, {
    rnd <- sample(edge_universe_ids, m, replace = FALSE)
    length(intersect(rnd, pathway_edge_ids))
  })
  
  (sum(hits >= obs) + 1) / (n_perm + 1)
}

run_enrichment_for_resource <- function(resource_module_edges, pathway_ref, edge_universe_ids,
                                        run_perm = FALSE, n_perm = 1000L) {
  module_edge_ids <- intersect(unique(resource_module_edges$edge_id), unique(edge_universe_ids))
  
  out <- pathway_ref %>%
    transmute(
      source_db, pathway_id, pathway_name,
      n_nodes_in_benchmark, n_pathway_edges_in_universe,
      pathway_edge_ids_in_universe
    ) %>%
    mutate(
      stat = map(pathway_edge_ids_in_universe,
                 ~ fisher_edge_enrichment(module_edge_ids, edge_universe_ids, .x))
    ) %>%
    unnest(stat) %>%
    mutate(
      p_adj = p.adjust(p_value, method = "BH"),
      neglog10_padj = -log10(pmax(p_adj, 1e-300)),
      edge_precision_in_pathway = overlap_edges / pmax(module_edges, 1),
      pathway_edge_coverage = overlap_edges / pmax(pathway_edges, 1),
      jaccard_edges = overlap_edges / pmax(module_edges + pathway_edges - overlap_edges, 1)
    )
  
  if (run_perm) {
    out <- out %>%
      mutate(
        empirical_p = map_dbl(pathway_edge_ids_in_universe,
                              ~ empirical_edge_p(module_edge_ids, edge_universe_ids, .x, n_perm = n_perm)),
        empirical_p_adj = p.adjust(empirical_p, method = "BH"),
        neglog10_empirical_padj = -log10(pmax(empirical_p_adj, 1e-300))
      )
  }
  
  out
}

# I/O for thresholded BioNet modules
read_module_edges_for_threshold <- function(thr, resource_key) {
  thr_candidates <- unique(c(
    formatC(thr, format = "f", digits = 2),
    formatC(thr, format = "f", digits = 1),
    as.character(thr)
  ))
  
  candidate_files <- file.path(
    "output",
    paste0("bionet_integrated_", resource_key),
    paste0("module_edges_fdr_", thr_candidates, ".csv")
  )
  
  hit_idx <- which(file.exists(candidate_files))
  if (length(hit_idx) == 0) {
    stop(glue(
      "Missing module edge file for resource='{resource_key}', threshold={thr}. ",
      "Tried: {paste(candidate_files, collapse = ' | ')}"
    ))
  }
  
  f <- candidate_files[hit_idx[1]]
  message(glue("  -> loading: {f}"))
  
  x <- readr::read_csv(f, show_col_types = FALSE)
  
  if (ncol(x) < 2) stop(glue("File has fewer than 2 columns (cannot parse edge list): {f}"))
  colnames(x)[1:2] <- c("geneA", "geneB")
  x
}

# LOAD / STANDARDIZE FULL INTERACTOMES (once)
caprinim_full <- canonicalize_edges_df(ppi_list$caprinim, "A", "B", "CAPRINIM")
biogrid_full  <- canonicalize_edges_df(ppi_list$biogrid,  "A", "B", "BioGRID")
string_full   <- canonicalize_edges_df(ppi_list$string,   "A", "B", "STRING")

full_edges_list <- list(CAPRINIM = caprinim_full, BioGRID = biogrid_full, STRING = string_full)
all_full_edges <- bind_rows(full_edges_list)

edge_universe_ids <- unique(all_full_edges$edge_id)
benchmark_genes <- unique(c(all_full_edges$gene1, all_full_edges$gene2))

message(glue("[FULL] CAPRINIM: {nrow(caprinim_full)} edges"))
message(glue("[FULL] BioGRID : {nrow(biogrid_full)} edges"))
message(glue("[FULL] STRING  : {nrow(string_full)} edges"))
message(glue("Benchmark edge universe (union of full interactomes): {length(edge_universe_ids)}"))
message(glue("Benchmark genes (from full interactomes): {length(benchmark_genes)}"))

caprinim_universe <- caprinim_full$edge_id
biogrid_universe  <- biogrid_full$edge_id
string_universe   <- string_full$edge_id

# LOAD ALL TOPOLOGY-CAPABLE PATHWAY RESOURCES (graphite), KEEP ONLY THOSE WITH EDGES
candidate_dbs <- discover_graphite_databases(verbose = TRUE)

combined_cache <- file.path("cache", paste0(
  "pathway_edges_topology_",
  graphite_species,
  "_ALL_AVAILABLE_DBS_collapse_", ifelse(collapse_pathway_direction, "T", "F"),
  "_cardio_", ifelse(restrict_to_cardio_pathways, "T", "F"),
  ".rds"
))

cache_meta_key <- list(
  species = graphite_species,
  candidate_dbs = sort(candidate_dbs),
  collapse_direction = collapse_pathway_direction,
  restrict_cardio = restrict_to_cardio_pathways
)

pathway_edges_topology <- NULL

if (file.exists(combined_cache)) {
  message(glue("Found combined cache: {combined_cache}"))
  obj <- readRDS(combined_cache)
  meta <- attr(obj, "meta")
  same <- !is.null(meta) &&
    identical(meta$species, cache_meta_key$species) &&
    identical(sort(meta$candidate_dbs), sort(cache_meta_key$candidate_dbs)) &&
    identical(meta$collapse_direction, cache_meta_key$collapse_direction) &&
    identical(meta$restrict_cardio, cache_meta_key$restrict_cardio)
  
  if (isTRUE(same)) {
    message("  -> cache metadata matches; using cached combined topology edges.")
    pathway_edges_topology <- obj
  } else {
    message("  -> cache metadata differs; recomputing combined topology edges.")
  }
}

if (is.null(pathway_edges_topology)) {
  db_reports <- vector("list", length(candidate_dbs))
  edges_list <- vector("list", length(candidate_dbs))
  
  for (i in seq_along(candidate_dbs)) {
    db <- candidate_dbs[i]
    cache_file <- file.path("cache", glue("pathway_edges_topology_{graphite_species}_{db}.rds"))
    
    message(glue("\n=== Trying pathway DB: {db} ({i}/{length(candidate_dbs)}) ==="))
    
    res <- tryCatch({
      ed <- load_graphite_topology_edges_db(
        species = graphite_species,
        pathway_db = db,
        collapse_direction = collapse_pathway_direction,
        restrict_cardio = restrict_to_cardio_pathways,
        cache_file = cache_file,
        use_cache = TRUE,
        verbose = TRUE
      )
      ed
    }, error = function(e) {
      attr(e, "db") <- db
      e
    })
    
    if (inherits(res, "error")) {
      db_reports[[i]] <- tibble(
        source_db = db,
        status = "FAILED",
        n_pathways = NA_integer_,
        n_edges = NA_integer_,
        error = conditionMessage(res)
      )
      edges_list[[i]] <- NULL
      message(glue("  !! FAILED {db}: {conditionMessage(res)}"))
      next
    }
    
    n_edges <- nrow(res)
    n_paths <- if (n_edges > 0) n_distinct(res$pathway_id) else 0L
    
    if (n_edges == 0) {
      db_reports[[i]] <- tibble(
        source_db = db,
        status = "EMPTY_EDGES",
        n_pathways = n_paths,
        n_edges = n_edges,
        error = NA_character_
      )
      edges_list[[i]] <- NULL
      message(glue("  -> {db}: returned 0 topology edges (skipping)."))
      next
    }
    
    # Keep only if there are edges and a reasonable number of pathways
    db_reports[[i]] <- tibble(
      source_db = db,
      status = "OK",
      n_pathways = n_paths,
      n_edges = n_edges,
      error = NA_character_
    )
    edges_list[[i]] <- res
    message(glue("  -> {db}: kept ({n_paths} pathways; {n_edges} pathway-edge entries)."))
  }
  
  db_report_tbl <- bind_rows(db_reports) %>%
    arrange(desc(status == "OK"), source_db)
  
  write_tsv(db_report_tbl, file.path(outdir, "tables", "pathway_db_load_report.tsv"))
  
  pathway_edges_topology <- bind_rows(compact(edges_list))
  if (nrow(pathway_edges_topology) == 0) {
    stop("No pathway databases yielded topology edges for this species/settings. Check graphite installation/species/cardio filter.")
  }
  
  attr(pathway_edges_topology, "meta") <- cache_meta_key
  saveRDS(pathway_edges_topology, combined_cache)
  message(glue("\nSaved combined cache: {combined_cache}"))
}

message(glue("\n=== Combined topology edges loaded: {nrow(pathway_edges_topology)} pathway-edge entries across ",
             "{n_distinct(pathway_edges_topology$source_db)} DB(s) and {n_distinct(pathway_edges_topology$pathway_id)} pathway IDs (DB-local) ==="))

write_tsv(
  pathway_edges_topology %>% select(source_db, pathway_id, pathway_name, source, target, edge_id),
  file.path(outdir, "tables", "pathway_topology_edges_raw.tsv")
)

pathway_ref <- build_pathway_topology_reference(
  pathway_edges_tbl = pathway_edges_topology,
  benchmark_genes = benchmark_genes,
  edge_universe_ids = edge_universe_ids,
  min_nodes = min_pathway_nodes_in_benchmark,
  min_edges = min_pathway_edges_in_universe
)

if (nrow(pathway_ref) == 0) {
  stop("No pathways overlap enough with your edge universe. Check species / gene IDs / mapping / cardio filter.")
}

write_tsv(
  pathway_ref %>%
    transmute(
      source_db,
      pathway_id,
      pathway_name = pretty_pathway_name(pathway_name),
      n_nodes_in_benchmark,
      n_pathway_edges_in_universe
    ) %>%
    arrange(source_db, desc(n_pathway_edges_in_universe)),
  file.path(outdir, "tables", "pathway_reference_summary_topology.tsv")
)

# CORE ANALYSIS LOOP
run_one_threshold <- function(thr) {
  thr_label <- fmt_thr(thr)
  message(glue("\n=== Running threshold {thr_label} ==="))
  
  caprinim_edges <- read_module_edges_for_threshold(thr, "caprinim")
  biogrid_edges  <- read_module_edges_for_threshold(thr, "biogrid")
  string_edges   <- read_module_edges_for_threshold(thr, "string")
  
  caprinim_module <- canonicalize_edges_df(caprinim_edges, "geneA", "geneB", "CAPRINIM")
  biogrid_module  <- canonicalize_edges_df(biogrid_edges,  "geneA", "geneB", "BioGRID")
  string_module   <- canonicalize_edges_df(string_edges,   "geneA", "geneB", "STRING")
  
  message(glue("[MODULE {thr_label}] CAPRINIM: {nrow(caprinim_module)} edges"))
  message(glue("[MODULE {thr_label}] BioGRID : {nrow(biogrid_module)} edges"))
  message(glue("[MODULE {thr_label}] STRING  : {nrow(string_module)} edges"))
  
  res_caprinim <- run_enrichment_for_resource(
    resource_module_edges = caprinim_module,
    pathway_ref = pathway_ref,
    edge_universe_ids = caprinim_universe,
    run_perm = run_empirical_permutations,
    n_perm = n_permutations
  ) %>% mutate(resource = "CAPRINIM")
  
  res_biogrid <- run_enrichment_for_resource(
    resource_module_edges = biogrid_module,
    pathway_ref = pathway_ref,
    edge_universe_ids = biogrid_universe,
    run_perm = run_empirical_permutations,
    n_perm = n_permutations
  ) %>% mutate(resource = "BioGRID")
  
  res_string <- run_enrichment_for_resource(
    resource_module_edges = string_module,
    pathway_ref = pathway_ref,
    edge_universe_ids = string_universe,
    run_perm = run_empirical_permutations,
    n_perm = n_permutations
  ) %>% mutate(resource = "STRING")
  
  all_results <- bind_rows(res_caprinim, res_biogrid, res_string) %>%
    mutate(
      threshold = thr,
      threshold_label = thr_label,
      pathway_label = paste0("[", toupper(source_db), "] ", pretty_pathway_name(pathway_name)),
      significant = p_adj < alpha
    )
  
  summary_resource <- all_results %>%
    group_by(threshold, threshold_label, resource) %>%
    summarise(
      n_pathways_tested = n(),
      n_significant = sum(significant, na.rm = TRUE),
      mean_neglog10_padj = mean(neglog10_padj[is.finite(neglog10_padj)], na.rm = TRUE),
      median_neglog10_padj = median(neglog10_padj[is.finite(neglog10_padj)], na.rm = TRUE),
      mean_jaccard = mean(jaccard_edges, na.rm = TRUE),
      mean_pathway_coverage = mean(pathway_edge_coverage, na.rm = TRUE),
      mean_edge_precision = mean(edge_precision_in_pathway, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      module_edges_in_universe = case_when(
        resource == "CAPRINIM" ~ length(intersect(caprinim_module$edge_id, edge_universe_ids)),
        resource == "BioGRID"  ~ length(intersect(biogrid_module$edge_id, edge_universe_ids)),
        resource == "STRING"   ~ length(intersect(string_module$edge_id, edge_universe_ids)),
        TRUE ~ NA_integer_
      )
    ) %>%
    arrange(threshold, desc(n_significant), desc(mean_neglog10_padj))
  
  wide_comp <- all_results %>%
    select(threshold, threshold_label, resource, source_db, pathway_id, pathway_label,
           p_adj, neglog10_padj, odds_ratio, overlap_edges, module_edges, pathway_edges,
           edge_precision_in_pathway, pathway_edge_coverage, jaccard_edges) %>%
    pivot_wider(
      names_from = resource,
      values_from = c(
        p_adj, neglog10_padj, odds_ratio, overlap_edges,
        edge_precision_in_pathway, pathway_edge_coverage, jaccard_edges
      )
    ) %>%
    mutate(
      delta_caprinim_vs_biogrid = neglog10_padj_CAPRINIM - neglog10_padj_BioGRID,
      delta_caprinim_vs_string  = neglog10_padj_CAPRINIM - neglog10_padj_STRING,
      mean_delta_caprinim = rowMeans(cbind(delta_caprinim_vs_biogrid, delta_caprinim_vs_string), na.rm = TRUE)
    ) %>%
    arrange(desc(mean_delta_caprinim))
  
  pathway_wins <- all_results %>%
    group_by(threshold, threshold_label, source_db, pathway_id, pathway_label) %>%
    mutate(rank_strength = min_rank(desc(neglog10_padj))) %>%
    filter(rank_strength == 1) %>%
    summarise(winners = paste(sort(unique(resource)), collapse = ","), .groups = "drop") %>%
    separate_rows(winners, sep = ",") %>%
    filter(winners != "") %>%
    count(threshold, threshold_label, winners, name = "n_pathway_wins") %>%
    rename(resource = winners) %>%
    arrange(threshold, desc(n_pathway_wins))
  
  cap_specific_hits <- all_results %>%
    select(threshold, threshold_label, resource, source_db, pathway_id, pathway_label, p_adj, neglog10_padj, overlap_edges) %>%
    pivot_wider(names_from = resource, values_from = c(p_adj, neglog10_padj, overlap_edges)) %>%
    mutate(
      cap_sig = p_adj_CAPRINIM < alpha,
      bg_sig  = p_adj_BioGRID < alpha,
      st_sig  = p_adj_STRING < alpha
    ) %>%
    filter(cap_sig, !bg_sig, !st_sig) %>%
    arrange(desc(neglog10_padj_CAPRINIM))
  
  list(
    all_results = all_results,
    summary_resource = summary_resource,
    wide_comp = wide_comp,
    pathway_wins = pathway_wins,
    cap_specific_hits = cap_specific_hits
  )
}

results_by_thr <- map(bionet_thresholds, run_one_threshold)

all_results_all <- map_dfr(results_by_thr, "all_results")
summary_resource_all <- map_dfr(results_by_thr, "summary_resource")
wide_comp_all <- map_dfr(results_by_thr, "wide_comp")
pathway_wins_all <- map_dfr(results_by_thr, "pathway_wins")
cap_specific_hits_all <- map_dfr(results_by_thr, "cap_specific_hits")

# T-TESTS / DISTRIBUTION COMPARISONS
run_pairwise_distribution_tests <- function(df_one_threshold) {
  wide <- df_one_threshold %>%
    select(source_db, pathway_id, pathway_label, resource, neglog10_padj) %>%
    pivot_wider(names_from = resource, values_from = neglog10_padj)
  
  paired_test_row <- function(x, y, label_x, label_y) {
    ok <- is.finite(x) & is.finite(y)
    x2 <- x[ok]
    y2 <- y[ok]
    
    if (length(x2) < 3) {
      return(tibble(
        group1 = label_x, group2 = label_y, n_pairs = length(x2),
        mean_group1 = mean(x2, na.rm = TRUE),
        mean_group2 = mean(y2, na.rm = TRUE),
        mean_diff = mean(x2 - y2, na.rm = TRUE),
        median_diff = median(x2 - y2, na.rm = TRUE),
        t_statistic = NA_real_, df = NA_real_, p_value = NA_real_,
        conf_low = NA_real_, conf_high = NA_real_,
        wilcox_p_value = NA_real_
      ))
    }
    
    tt <- t.test(x2, y2, paired = TRUE, alternative = "two.sided")
    wt <- suppressWarnings(wilcox.test(x2, y2, paired = TRUE, exact = FALSE))
    
    tibble(
      group1 = label_x,
      group2 = label_y,
      n_pairs = length(x2),
      mean_group1 = mean(x2, na.rm = TRUE),
      mean_group2 = mean(y2, na.rm = TRUE),
      mean_diff = mean(x2 - y2, na.rm = TRUE),
      median_diff = median(x2 - y2, na.rm = TRUE),
      t_statistic = unname(tt$statistic),
      df = unname(tt$parameter),
      p_value = tt$p.value,
      conf_low = tt$conf.int[1],
      conf_high = tt$conf.int[2],
      wilcox_p_value = wt$p.value
    )
  }
  
  bind_rows(
    paired_test_row(wide$CAPRINIM, wide$BioGRID, "CAPRINIM", "BioGRID"),
    paired_test_row(wide$CAPRINIM, wide$STRING,  "CAPRINIM", "STRING")
  )
}

ttest_by_threshold <- all_results_all %>%
  group_by(threshold, threshold_label) %>%
  group_modify(~ run_pairwise_distribution_tests(.x)) %>%
  ungroup() %>%
  mutate(
    p_adj_bh = p.adjust(p_value, method = "BH"),
    wilcox_p_adj_bh = p.adjust(wilcox_p_value, method = "BH")
  )

paired_diffs_all <- all_results_all %>%
  select(threshold, threshold_label, source_db, pathway_id, pathway_label, resource, neglog10_padj) %>%
  pivot_wider(names_from = resource, values_from = neglog10_padj) %>%
  mutate(
    diff_cap_vs_biogrid = CAPRINIM - BioGRID,
    diff_cap_vs_string  = CAPRINIM - STRING
  ) %>%
  select(threshold, threshold_label, source_db, pathway_id, pathway_label, diff_cap_vs_biogrid, diff_cap_vs_string) %>%
  pivot_longer(
    cols = c(diff_cap_vs_biogrid, diff_cap_vs_string),
    names_to = "comparison",
    values_to = "delta_neglog10_padj"
  ) %>%
  mutate(
    comparison = recode(comparison,
                        diff_cap_vs_biogrid = "CAPRINIM - BioGRID",
                        diff_cap_vs_string  = "CAPRINIM - STRING")
  )

# WRITE TABLES
write_tsv(all_results_all, file.path(outdir, "tables", "edge_topology_enrichment_all_results_threshold_0p05.tsv"))
write_tsv(summary_resource_all, file.path(outdir, "tables", "comparison_summary_by_resource_threshold_0p05.tsv"))
write_tsv(wide_comp_all, file.path(outdir, "tables", "comparison_pathway_wide_topology_threshold_0p05.tsv"))
write_tsv(pathway_wins_all, file.path(outdir, "tables", "comparison_pathway_wins_threshold_0p05.tsv"))
write_tsv(cap_specific_hits_all, file.path(outdir, "tables", "caprinim_specific_topology_hits_threshold_0p05.tsv"))
write_tsv(ttest_by_threshold, file.path(outdir, "tables", "ttest_neglog10padj_caprinim_vs_baselines_threshold_0p05.tsv"))
write_tsv(paired_diffs_all, file.path(outdir, "tables", "paired_differences_neglog10padj_threshold_0p05.tsv"))

ttest_summary_pretty <- ttest_by_threshold %>%
  mutate(
    threshold_label = paste0("FDR ", threshold_label),
    comparison = paste(group1, "vs", group2),
    mean_diff = round(mean_diff, 3),
    t_statistic = round(t_statistic, 3),
    p_value = signif(p_value, 3),
    p_adj_bh = signif(p_adj_bh, 3),
    wilcox_p_value = signif(wilcox_p_value, 3),
    wilcox_p_adj_bh = signif(wilcox_p_adj_bh, 3)
  )

# PLOTS
resource_pal <- c(
  "CAPRINIM" = "#7C3AED",
  "BioGRID"  = "#0284C7",
  "STRING"   = "#EA580C"
)

# UPDATED: make paired-delta comparison colors match your requested BioGRID/STRING colors
comparison_pal <- c(
  "CAPRINIM - BioGRID" = "#0284C7",
  "CAPRINIM - STRING"  = "#EA580C"
)

thr_pal <- c("0.05" = "#BE123C")

base_theme <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "#E5E7EB"),
    panel.grid.major.y = element_line(color = "#EEF2FF"),
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(color = "grey20"),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )

# COMBINED PLOTS
p_counts <- if (nrow(summary_resource_all) == 0) {
  empty_plot("Significant pathway enrichments", "summary_resource_all is empty")
} else {
  summary_resource_all %>%
    mutate(threshold_label = factor(threshold_label, levels = fmt_thr(bionet_thresholds))) %>%
    ggplot(aes(x = threshold_label, y = n_significant, fill = resource)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.65, color = "white", linewidth = 0.25) +
    geom_text(aes(label = n_significant),
              position = position_dodge(width = 0.75), vjust = -0.25, size = 3.4, fontface = "bold") +
    scale_fill_manual(values = resource_pal, drop = FALSE) +
    labs(
      title = "Significant pathway enrichments at BioNet threshold 0.05",
      subtitle = glue("Edge-level topology enrichment (BH FDR < {alpha}); cardio-focused filter; all graphite DBs with topology edges"),
      x = "BioNet module edge FDR threshold",
      y = "Number of significant pathways",
      fill = "Resource"
    ) +
    expand_limits(y = max(summary_resource_all$n_significant, na.rm = TRUE) * 1.12) +
    base_theme
}

p_mean_enrich <- if (nrow(summary_resource_all) == 0) {
  empty_plot("Mean enrichment score", "summary_resource_all is empty")
} else {
  summary_resource_all %>%
    ggplot(aes(x = resource, y = mean_neglog10_padj, fill = resource)) +
    geom_col(width = 0.7, color = "white", linewidth = 0.25) +
    geom_text(aes(label = sprintf("%.2f", mean_neglog10_padj)), vjust = -0.3, size = 3.3, fontface = "bold") +
    scale_fill_manual(values = resource_pal, drop = FALSE) +
    labs(
      title = "Mean enrichment score at threshold 0.05",
      subtitle = "Metric: mean -log10(BH FDR) across matched pathways (all included DBs combined)",
      x = NULL,
      y = "Mean -log10(BH FDR)"
    ) +
    base_theme +
    theme(legend.position = "none")
}

dist_annot_df <- if (nrow(ttest_by_threshold) == 0) {
  tibble()
} else {
  ttest_by_threshold %>%
    mutate(
      threshold_label = factor(threshold_label, levels = fmt_thr(bionet_thresholds)),
      comparison_lab = paste(group1, "vs", group2),
      stars = sig_stars(p_adj_bh),
      label = glue("{comparison_lab}: t={round(t_statistic, 2)}, p={signif(p_value, 3)}, BH={signif(p_adj_bh, 3)} ({stars})")
    ) %>%
    group_by(threshold_label) %>%
    mutate(line_id = row_number()) %>%
    ungroup()
}

p_dist <- if (nrow(all_results_all) == 0) {
  empty_plot("Distribution of enrichment scores", "all_results_all is empty")
} else {
  p_base <- all_results_all %>%
    mutate(threshold_label = factor(threshold_label, levels = fmt_thr(bionet_thresholds))) %>%
    ggplot(aes(x = resource, y = neglog10_padj, fill = resource, color = resource)) +
    geom_violin(alpha = 0.22, trim = FALSE, linewidth = 0.7) +
    geom_boxplot(width = 0.17, outlier.shape = NA, alpha = 0.95, linewidth = 0.45, color = "#1F2937", fill = "white") +
    geom_jitter(width = 0.12, alpha = 0.16, size = 0.75, show.legend = FALSE) +
    facet_wrap(~ threshold_label, nrow = 1, drop = FALSE) +
    scale_fill_manual(values = resource_pal, drop = FALSE) +
    scale_color_manual(values = resource_pal, drop = FALSE) +
    labs(
      title = "Distribution of enrichment scores by resource (threshold 0.05)",
      subtitle = "Each point is a pathway (across all included topology-capable graphite DBs)",
      x = NULL,
      y = "-log10(BH FDR)"
    ) +
    base_theme +
    theme(axis.text.x = element_text(face = "bold"))
  
  if (nrow(dist_annot_df) > 0) {
    y_max_by_thr <- all_results_all %>%
      group_by(threshold_label) %>%
      summarise(ymax = max(neglog10_padj[is.finite(neglog10_padj)], na.rm = TRUE), .groups = "drop") %>%
      mutate(ymax = ifelse(is.finite(ymax), ymax, 1))
    
    ann <- dist_annot_df %>%
      left_join(y_max_by_thr, by = "threshold_label") %>%
      mutate(x = 2, y = ymax * (1.07 + 0.08 * (line_id - 1)))
    
    p_base +
      geom_text(
        data = ann,
        aes(x = x, y = y, label = label),
        inherit.aes = FALSE,
        size = 3.2,
        fontface = "bold",
        color = "#111827"
      ) +
      expand_limits(y = max(ann$y, na.rm = TRUE) * 1.05)
  } else {
    p_base
  }
}

p_delta_dist <- if (nrow(paired_diffs_all) == 0 || all(!is.finite(paired_diffs_all$delta_neglog10_padj))) {
  empty_plot("Paired pathway-level differences", "No finite paired differences available")
} else {
  paired_diffs_all %>%
    mutate(threshold_label = factor(threshold_label, levels = fmt_thr(bionet_thresholds))) %>%
    ggplot(aes(x = comparison, y = delta_neglog10_padj, fill = comparison)) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey40", linewidth = 0.6) +
    geom_violin(alpha = 0.28, trim = FALSE, linewidth = 0.7, color = NA) +
    geom_boxplot(width = 0.16, outlier.shape = NA, fill = "white", color = "#1F2937", linewidth = 0.45) +
    geom_jitter(aes(color = comparison), width = 0.10, alpha = 0.16, size = 0.75, show.legend = FALSE) +
    facet_wrap(~ threshold_label, nrow = 1, drop = FALSE) +
    scale_fill_manual(values = comparison_pal, drop = FALSE) +
    scale_color_manual(values = comparison_pal, drop = FALSE) +
    labs(
      title = "Paired pathway-level differences in enrichment scores",
      subtitle = "Positive values indicate stronger CAPRINIM enrichment than the comparator",
      x = NULL,
      y = expression(Delta * " -log10(BH FDR)")
    ) +
    base_theme +
    theme(axis.text.x = element_text(face = "bold"))
}

scatter_bg_all <- all_results_all %>%
  select(threshold, threshold_label, source_db, pathway_id, pathway_label, resource, neglog10_padj) %>%
  pivot_wider(names_from = resource, values_from = neglog10_padj) %>%
  filter(is.finite(CAPRINIM), is.finite(BioGRID))

p_scatter_bg <- if (nrow(scatter_bg_all) == 0) {
  empty_plot("CAPRINIM vs BioGRID enrichment (matched pathways)", "No matched finite pathway scores")
} else {
  ggplot(scatter_bg_all, aes(x = BioGRID, y = CAPRINIM)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey45") +
    geom_point(color = resource_pal["CAPRINIM"], alpha = 0.42, size = 1.6) +
    facet_wrap(~ threshold_label, nrow = 1, drop = FALSE) +
    labs(
      title = "CAPRINIM vs BioGRID enrichment (matched pathways)",
      subtitle = "Points above diagonal favor CAPRINIM",
      x = "BioGRID: -log10(BH FDR)",
      y = "CAPRINIM: -log10(BH FDR)"
    ) +
    coord_equal() +
    base_theme
}

scatter_st_all <- all_results_all %>%
  select(threshold, threshold_label, source_db, pathway_id, pathway_label, resource, neglog10_padj) %>%
  pivot_wider(names_from = resource, values_from = neglog10_padj) %>%
  filter(is.finite(CAPRINIM), is.finite(STRING))

p_scatter_st <- if (nrow(scatter_st_all) == 0) {
  empty_plot("CAPRINIM vs STRING enrichment (matched pathways)", "No matched finite pathway scores")
} else {
  ggplot(scatter_st_all, aes(x = STRING, y = CAPRINIM)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey45") +
    geom_point(color = resource_pal["STRING"], alpha = 0.42, size = 1.6) +
    facet_wrap(~ threshold_label, nrow = 1, drop = FALSE) +
    labs(
      title = "CAPRINIM vs STRING enrichment (matched pathways)",
      subtitle = "Points above diagonal favor CAPRINIM",
      x = "STRING: -log10(BH FDR)",
      y = "CAPRINIM: -log10(BH FDR)"
    ) +
    coord_equal() +
    base_theme
}

top_delta_combined <- wide_comp_all %>%
  mutate(threshold_label = factor(threshold_label, levels = fmt_thr(bionet_thresholds))) %>%
  filter(is.finite(mean_delta_caprinim)) %>%
  group_by(threshold, threshold_label) %>%
  slice_max(order_by = mean_delta_caprinim, n = 12, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(pathway_label = str_trunc(pathway_label, width = 60))

p_top_delta <- if (nrow(top_delta_combined) == 0) {
  empty_plot("Top pathways where CAPRINIM outperforms BioGRID/STRING",
             "No finite mean delta values were available")
} else {
  top_delta_combined <- top_delta_combined %>%
    mutate(pathway_label_f = paste0(pathway_label, "  [", threshold_label, "]")) %>%
    group_by(threshold_label) %>%
    mutate(pathway_label_f = fct_reorder(pathway_label_f, mean_delta_caprinim)) %>%
    ungroup()
  
  ggplot(top_delta_combined, aes(x = mean_delta_caprinim, y = pathway_label_f, fill = threshold_label)) +
    geom_col(width = 0.76, color = "white", linewidth = 0.25) +
    facet_wrap(~ threshold_label, scales = "free_y", ncol = 1, drop = FALSE) +
    scale_fill_manual(values = thr_pal, guide = "none", drop = FALSE) +
    labs(
      title = "Top pathways where CAPRINIM outperforms BioGRID/STRING",
      subtitle = "Ranking by mean delta = avg(CAPRINIM-BioGRID, CAPRINIM-STRING)",
      x = "Mean delta -log10(BH FDR)",
      y = NULL
    ) +
    base_theme +
    theme(axis.text.y = element_text(size = 8))
}

tt_plot_df <- ttest_by_threshold %>%
  mutate(
    comparison = factor(paste(group1, "vs", group2),
                        levels = c("CAPRINIM vs BioGRID", "CAPRINIM vs STRING")),
    threshold_label = factor(threshold_label, levels = fmt_thr(bionet_thresholds))
  ) %>%
  filter(is.finite(mean_diff), is.finite(conf_low), is.finite(conf_high))

p_ttest <- if (nrow(tt_plot_df) == 0) {
  empty_plot("Paired t-test effect sizes", "No valid paired comparisons were available")
} else {
  ggplot(tt_plot_df, aes(x = threshold_label, y = mean_diff, color = comparison, group = comparison)) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey45") +
    geom_errorbar(aes(ymin = conf_low, ymax = conf_high),
                  width = 0.08, linewidth = 0.9,
                  position = position_dodge(width = 0.25)) +
    geom_point(size = 3.0, position = position_dodge(width = 0.25)) +
    geom_line(linewidth = 1.05, position = position_dodge(width = 0.25)) +
    scale_color_manual(values = comparison_pal, drop = FALSE) +
    labs(
      title = "Paired t-test effect sizes",
      subtitle = "Mean paired difference in pathway enrichment score (CAPRINIM - comparator) with 95% CI",
      x = "BioNet module edge FDR threshold",
      y = "Mean paired difference in -log10(BH FDR)",
      color = "Comparison"
    ) +
    base_theme
}

wins_plot_df <- pathway_wins_all %>%
  mutate(
    threshold_label = factor(threshold_label, levels = fmt_thr(bionet_thresholds)),
    resource = factor(resource, levels = c("CAPRINIM", "BioGRID", "STRING"))
  )

p_wins <- if (nrow(wins_plot_df) == 0) {
  empty_plot("Number of pathway wins by resource", "No pathway wins available")
} else {
  ggplot(wins_plot_df, aes(x = resource, y = n_pathway_wins, fill = resource)) +
    geom_col(width = 0.68, color = "white", linewidth = 0.25) +
    geom_text(aes(label = n_pathway_wins), vjust = -0.25, size = 3.4, fontface = "bold") +
    facet_wrap(~ threshold_label, nrow = 1, drop = FALSE) +
    scale_fill_manual(values = resource_pal, drop = FALSE) +
    labs(
      title = "Number of pathway wins by resource",
      subtitle = "A win = strongest enrichment (-log10 BH FDR); ties counted for all winners",
      x = NULL,
      y = "Number of wins"
    ) +
    base_theme +
    theme(legend.position = "none")
}

# COMBINED FIGURES
fig_main <- (p_counts | p_mean_enrich) / (p_dist / p_delta_dist) +
  plot_annotation(
    title = "Edge-topology enrichment benchmark (BioNet threshold 0.05)",
    subtitle = "CAPRINIM vs BioGRID vs STRING; all topology-capable graphite pathway DBs (filtered)",
    theme = theme(
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 11, color = "grey25")
    )
  )

fig_scatter <- p_scatter_bg / p_scatter_st
fig_cap_adv <- p_top_delta / p_ttest + plot_layout(heights = c(2.2, 1))
fig_wins <- p_wins

cat("\n--- Plot diagnostics ---\n")
cat("Rows in all_results_all: ", nrow(all_results_all), "\n")
cat("Rows in summary_resource_all: ", nrow(summary_resource_all), "\n")
cat("Rows in wide_comp_all: ", nrow(wide_comp_all), "\n")
cat("Rows in top_delta_combined: ", nrow(top_delta_combined), "\n")
cat("Rows in ttest_by_threshold: ", nrow(ttest_by_threshold), "\n")
cat("Rows in tt_plot_df (finite): ", nrow(tt_plot_df), "\n")
cat("Rows in scatter_bg_all: ", nrow(scatter_bg_all), "\n")
cat("Rows in scatter_st_all: ", nrow(scatter_st_all), "\n")
cat("Rows in pathway_wins_all: ", nrow(pathway_wins_all), "\n")

# HEATMAP OF TOP PATHWAYS (clustered)
top_n_per_threshold <- 25

top_paths_multi <- all_results_all %>%
  group_by(threshold, threshold_label, source_db, pathway_id, pathway_label) %>%
  summarise(max_enrich = max(neglog10_padj, na.rm = TRUE), .groups = "drop") %>%
  mutate(max_enrich = ifelse(is.infinite(max_enrich), NA_real_, max_enrich)) %>%
  filter(is.finite(max_enrich)) %>%
  group_by(threshold, threshold_label) %>%
  slice_max(order_by = max_enrich, n = top_n_per_threshold, with_ties = FALSE) %>%
  ungroup()

p_heatmap_multi <- if (nrow(top_paths_multi) == 0) {
  empty_plot(
    title = "Top pathway enrichments (clustered heatmap)",
    subtitle = "No finite enrichment values available to build heatmap"
  )
} else {
  heatmap_df_multi_raw <- top_paths_multi %>%
    select(threshold, threshold_label, source_db, pathway_id, pathway_label) %>%
    inner_join(
      all_results_all %>%
        select(threshold, threshold_label, resource, source_db, pathway_id, pathway_label, neglog10_padj),
      by = c("threshold", "threshold_label", "source_db", "pathway_id", "pathway_label")
    ) %>%
    mutate(
      pathway_label = str_trunc(pathway_label, 45),
      threshold_label = factor(threshold_label, levels = fmt_thr(bionet_thresholds))
    )
  
  clustered_heatmap_list <- heatmap_df_multi_raw %>% group_split(threshold_label)
  
  clustered_heatmap_df <- map_dfr(clustered_heatmap_list, function(df_thr) {
    ord <- cluster_heatmap_order(df_thr, row_col = "pathway_label", col_col = "resource", value_col = "neglog10_padj")
    df_thr %>%
      mutate(
        pathway_label = factor(pathway_label, levels = rev(unique(ord$row_order))),
        resource = factor(resource, levels = ord$col_order)
      )
  })
  
  ggplot(clustered_heatmap_df, aes(x = resource, y = pathway_label, fill = neglog10_padj)) +
    geom_tile(color = "white", linewidth = 0.2) +
    facet_wrap(~ threshold_label, nrow = 1, scales = "free_y", drop = FALSE) +
    scale_fill_gradientn(colors = c("#F8FAFC", "#BFDBFE", "#60A5FA", "#2563EB", "#1D4ED8", "#1E1B4B")) +
    labs(
      title = "Top pathway enrichments across threshold(s) (clustered heatmap)",
      subtitle = glue("Top {top_n_per_threshold} pathways per threshold by maximum enrichment across resources; rows clustered"),
      x = NULL, y = NULL, fill = "-log10(BH FDR)"
    ) +
    base_theme +
    theme(axis.text.y = element_text(size = 7))
}

# SAVE FIGURES AS PDF
safe_ggsave(file.path(outdir, "plots", "combined_main_summary_distributions.pdf"),
            fig_main, width = 16, height = 14, dpi = 320)

safe_ggsave(file.path(outdir, "plots", "combined_scatter_comparisons.pdf"),
            fig_scatter, width = 14, height = 10, dpi = 320)

safe_ggsave(file.path(outdir, "plots", "combined_caprinim_advantage_and_ttests.pdf"),
            fig_cap_adv, width = 14, height = 13, dpi = 320)

safe_ggsave(file.path(outdir, "plots", "plot_significant_counts_threshold_0p05.pdf"),
            p_counts, width = 10, height = 5.5, dpi = 320)

safe_ggsave(file.path(outdir, "plots", "plot_distribution_enrichment_scores_by_threshold.pdf"),
            p_dist, width = 16, height = 6.6, dpi = 320)

safe_ggsave(file.path(outdir, "plots", "plot_paired_deltas_caprinim_vs_baselines.pdf"),
            p_delta_dist, width = 12, height = 5.8, dpi = 320)

safe_ggsave(file.path(outdir, "plots", "plot_ttest_effectsizes_threshold_0p05.pdf"),
            p_ttest, width = 10, height = 5.5, dpi = 320)

safe_ggsave(file.path(outdir, "plots", "plot_pathway_wins_by_resource_threshold_0p05.pdf"),
            p_wins, width = 10, height = 5.8, dpi = 320)

safe_ggsave(file.path(outdir, "plots", "combined_heatmap_top_pathways_across_thresholds.pdf"),
            p_heatmap_multi, width = 9, height = 16, dpi = 320)

# CONSOLE SUMMARY
message("\n=== DONE: edge-level topology enrichment benchmark (threshold 0.05) ===")
message(glue("Thresholds screened: {paste(fmt_thr(bionet_thresholds), collapse = ', ')}"))
message(glue("Output directory: {normalizePath(outdir)}"))

cat("\n--- Resource summary ---\n")
print(summary_resource_all)

cat("\n--- Paired t-tests on enrichment scores (-log10 BH FDR) ---\n")
print(ttest_summary_pretty)

top_caprinim_by_thr <- all_results_all %>%
  filter(resource == "CAPRINIM") %>%
  arrange(threshold, p_adj, desc(overlap_edges)) %>%
  group_by(threshold_label) %>%
  slice_head(n = 15) %>%
  ungroup() %>%
  select(threshold_label, source_db, pathway_label, p_adj, neglog10_padj, odds_ratio, overlap_edges, pathway_edges, pathway_edge_coverage)

cat("\n--- Top CAPRINIM pathways (preview) ---\n")
print(top_caprinim_by_thr)

message(glue("\nPathway DB load report saved to: {file.path(outdir, 'tables', 'pathway_db_load_report.tsv')}"))