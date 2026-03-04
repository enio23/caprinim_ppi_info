#! /usr/bin/env Rscript

set.seed(1234)

library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(igraph)
library(BioNet)

dir.create("output", showWarnings = FALSE)

# Simple progress / logging helpers
ts_now <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

log_msg <- function(..., logfile = NULL) {
  msg <- paste0("[", ts_now(), "] ", paste(..., collapse = ""))
  message(msg)
  if (!is.null(logfile)) {
    cat(msg, "\n", file = logfile, append = TRUE)
  }
}

# Load data
logfile_global <- file.path("output", "bionet_runlog_integrated.txt")
log_msg("Starting BioNet *integrated* pipeline", logfile = logfile_global)

load(file = "/beegfs/prj/Dewenter_TAC_Backs_lab/Analyses_Gjerga_Questions/all_results.RData")
load(file = "output/ppi_list.RData")
names(ppi_list) <- c("caprinim", "biogrid", "string")

log_msg("Loaded all_results and ppi_list with PPIs: ", paste(names(ppi_list), collapse = ", "),
        logfile = logfile_global)

# Build df exactly as you do (common genes between proteomics and RNA)
dpa <- all_results$dpa_expr[which(all_results$dpa_expr$tp.comparison == "TAC (Combined) vs Sham (Combined)"), ]

dge <- readRDS(file = "/beegfs/prj/Dewenter_TAC_Backs_lab/Analyses_Gjerga/sequencing_analysis_cdna_assembly/with_ctrl_0/output/QlfList.rds")
dge <- dge$Condition

common_genes <- intersect(x = dpa$SYMBOL, y = dge$SYMBOL)
log_msg("Common genes (intersection proteomics x RNA): ", length(common_genes), logfile = logfile_global)

df <- matrix(nrow = length(common_genes), ncol = 5)
for(ii in seq_along(common_genes)){
  df[ii, 1] <- common_genes[ii]
  df[ii, 2] <- dpa$logFC[which(dpa$SYMBOL == common_genes[ii])[1]]
  df[ii, 3] <- dpa$p.val[which(dpa$SYMBOL == common_genes[ii])[1]] 
  df[ii, 4] <- dge$logFC[which(dge$SYMBOL == common_genes[ii])[1]]
  df[ii, 5] <- dge$PValue[which(dge$SYMBOL == common_genes[ii])[1]]
}
colnames(df) <- c("symbol", "prot_logfc", "prot_padj", "gex_logfc", "gex_padj")
df <- as.data.frame(df)

df$prot_logfc <- as.numeric(df$prot_logfc)
df$prot_padj  <- as.numeric(df$prot_padj)
df$gex_logfc  <- as.numeric(df$gex_logfc)
df$gex_padj   <- as.numeric(df$gex_padj)

# Clean & prep omics table
omics <- df %>%
  transmute(
    symbol = as.character(symbol),
    prot_logfc = as.numeric(prot_logfc),
    prot_padj  = as.numeric(prot_padj),
    gex_logfc  = as.numeric(gex_logfc),
    gex_padj   = as.numeric(gex_padj)
  ) %>%
  filter(!is.na(symbol), symbol != "") %>%
  distinct(symbol, .keep_all = TRUE)

# normalize p-values: invalid -> 1
omics$prot_padj[is.na(omics$prot_padj) | omics$prot_padj <= 0 | omics$prot_padj > 1] <- 1
omics$gex_padj [is.na(omics$gex_padj)  | omics$gex_padj  <= 0 | omics$gex_padj  > 1] <- 1

log_msg("Omics table rows: ", nrow(omics),
        " | prot_p<0.05: ", sum(omics$prot_padj < 0.05),
        " | gex_p<0.05: ", sum(omics$gex_padj < 0.05),
        logfile = logfile_global)

# Integrate proteomics + RNA evidence using BioNet::aggrPvals
pvals_mat <- cbind(prot = omics$prot_padj, gex = omics$gex_padj)
rownames(pvals_mat) <- omics$symbol

pval_integrated <- BioNet::aggrPvals(pvals_mat, order = 2, plot = FALSE)

omics$int_pval <- as.numeric(pval_integrated[omics$symbol])

log_msg("Integrated p-values computed with aggrPvals(order=2).",
        " | int_p<0.05: ", sum(omics$int_pval < 0.05, na.rm = TRUE),
        logfile = logfile_global)

# run BioNet for ONE PPI using integrated p-values
run_bionet_integrated_ppi <- function(ppi_df,
                                      ppi_name,
                                      omics,
                                      fdr_vec = c(0.05),
                                      keep_lcc = TRUE,
                                      restrict_to_measured = TRUE,
                                      measured_mode = c("any","both"),
                                      logfile = NULL) {
  
  measured_mode <- match.arg(measured_mode)
  
  out_prefix <- "bionet_integrated"
  log_msg("---- ", out_prefix, " | PPI=", ppi_name, " | start ----", logfile = logfile)
  
  # Build graph from A/B columns
  edges <- ppi_df %>%
    dplyr::transmute(from = as.character(A), to = as.character(B)) %>%
    dplyr::filter(!is.na(from), !is.na(to), from != "", to != "", from != to)
  
  g <- igraph::graph_from_data_frame(edges, directed = FALSE)
  g <- igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
  
  if (igraph::vcount(g) == 0) {
    log_msg("PPI '", ppi_name, "' produced an empty graph. Skipping.", logfile = logfile)
    return(invisible(FALSE))
  }
  
  log_msg("Graph built: nodes=", igraph::vcount(g), " edges=", igraph::ecount(g), logfile = logfile)
  
  if (restrict_to_measured) {
    if (measured_mode == "any") {
      measured <- omics$symbol[(omics$prot_padj < 1) | (omics$gex_padj < 1)]
    } else {
      measured <- omics$symbol[(omics$prot_padj < 1) & (omics$gex_padj < 1)]
    }
    measured <- intersect(igraph::V(g)$name, measured)
    g <- igraph::induced_subgraph(g, vids = measured)
    g <- igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
    
    log_msg("After restrict_to_measured (", measured_mode, "): nodes=", igraph::vcount(g),
            " edges=", igraph::ecount(g), logfile = logfile)
  }
  
  # Keep largest connected component (as in tutorial: focus on connected modules) :contentReference[oaicite:7]{index=7}
  if (keep_lcc && igraph::vcount(g) > 0) {
    comp <- igraph::components(g)
    g <- igraph::induced_subgraph(g, which(comp$membership == which.max(comp$csize)))
    log_msg("After LCC: nodes=", igraph::vcount(g),
            " edges=", igraph::ecount(g), logfile = logfile)
  }
  
  if (igraph::vcount(g) == 0) {
    log_msg("Graph has 0 nodes after filtering. Skipping.", logfile = logfile)
    return(invisible(FALSE))
  }
  
  out_dir <- file.path("output", paste0(out_prefix, "_", ppi_name))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Map integrated p-values onto nodes (default 1)
  nodes <- igraph::V(g)$name
  pvals <- rep(1, length(nodes))
  names(pvals) <- nodes
  
  om_sub <- omics %>% dplyr::filter(symbol %in% nodes)
  if (nrow(om_sub) > 0) {
    pvals[om_sub$symbol] <- om_sub$int_pval
  }
  
  log_msg("Nodes with omics evidence in this graph: ", nrow(om_sub), "/", length(nodes),
          logfile = logfile)
  
  # Fit BUM on integrated p-value distribution :contentReference[oaicite:8]{index=8} :contentReference[oaicite:9]{index=9}
  fb <- suppressWarnings(BioNet::fitBumModel(pvals, plot = FALSE))
  
  # For plotting, BioNet::plotModule only takes one diff.expr vector.
  # Here we choose RNA logFC (gex_logfc) to color nodes, but we also export both logFCs in CSV.
  diff_expr <- rep(0, length(nodes))
  names(diff_expr) <- nodes
  if (nrow(om_sub) > 0) diff_expr[om_sub$symbol] <- om_sub$gex_logfc
  
  # Save meta
  meta <- list(
    ppi_name = ppi_name,
    integration = list(
      method = "BioNet::aggrPvals",
      order = 2,
      inputs = c("prot_padj", "gex_padj")
    ),
    nodes_total = igraph::vcount(g),
    edges_total = igraph::ecount(g),
    keep_lcc = keep_lcc,
    restrict_to_measured = restrict_to_measured,
    measured_mode = measured_mode
  )
  saveRDS(meta, file = file.path(out_dir, "meta.rds"))
  
  # Run per FDR: score nodes, run heuristic Heinz, export module :contentReference[oaicite:10]{index=10} :contentReference[oaicite:11]{index=11}
  for (fdr in fdr_vec) {
    tag <- paste0("fdr_", format(fdr, scientific = FALSE))
    log_msg("  [", out_prefix, " | ", ppi_name, "] scoring @ FDR=", fdr, logfile = logfile)
    
    scores <- BioNet::scoreNodes(network = g, fb = fb, fdr = fdr)
    saveRDS(scores, file = file.path(out_dir, paste0("scores_", tag, ".rds")))
    
    npos <- sum(scores > 0, na.rm = TRUE)
    log_msg("  [", out_prefix, " | ", ppi_name, "] positive nodes: ", npos, "/", length(scores),
            logfile = logfile)
    
    writeLines(
      sprintf("ppi=%s  integrated=aggrPvals(order=2)  fdr=%s  positive_nodes=%d  total_nodes=%d",
              ppi_name, fdr, npos, length(scores)),
      con = file.path(out_dir, paste0("score_summary_", tag, ".txt"))
    )
    
    if (npos == 0) {
      empty_nodes <- data.frame(symbol = character(0), stringsAsFactors = FALSE) %>%
        dplyr::left_join(omics, by = c("symbol" = "symbol"))
      utils::write.csv(empty_nodes,
                       file = file.path(out_dir, paste0("module_nodes_", tag, ".csv")),
                       row.names = FALSE)
      
      empty_edges <- data.frame(from = character(0), to = character(0), stringsAsFactors = FALSE)
      utils::write.csv(empty_edges,
                       file = file.path(out_dir, paste0("module_edges_", tag, ".csv")),
                       row.names = FALSE)
      
      saveRDS(igraph::make_empty_graph(), file = file.path(out_dir, paste0("module_", tag, ".rds")))
      
      writeLines(
        sprintf("No positive nodes at fdr=%s. Empty module.", fdr),
        con = file.path(out_dir, paste0("NOTE_", tag, ".txt"))
      )
      next
    }
    
    log_msg("  [", out_prefix, " | ", ppi_name, "] running runFastHeinz()", logfile = logfile)
    module <- BioNet::runFastHeinz(g, scores)
    saveRDS(module, file = file.path(out_dir, paste0("module_", tag, ".rds")))
    
    # Export nodes/edges (include integrated pval + both logFCs)
    mod_nodes <- data.frame(symbol = igraph::V(module)$name, stringsAsFactors = FALSE) %>%
      dplyr::left_join(omics, by = "symbol")
    utils::write.csv(mod_nodes,
                     file = file.path(out_dir, paste0("module_nodes_", tag, ".csv")),
                     row.names = FALSE)
    
    mod_edges <- igraph::as_data_frame(module, what = "edges")
    utils::write.csv(mod_edges,
                     file = file.path(out_dir, paste0("module_edges_", tag, ".csv")),
                     row.names = FALSE)
    
    # Plot (colored by gex_logfc here; change to prot_logfc if you prefer)
    log_msg("  [", out_prefix, " | ", ppi_name, "] writing PDF plot", logfile = logfile)
    try({
      grDevices::pdf(file = file.path(out_dir, paste0("module_plot_", tag, ".pdf")),
                     width = 9, height = 7)
      BioNet::plotModule(module, scores = scores, diff.expr = diff_expr)
      grDevices::dev.off()
    }, silent = TRUE)
    
    # GraphML export
    log_msg("  [", out_prefix, " | ", ppi_name, "] writing GraphML", logfile = logfile)
    igraph::write_graph(module,
                        file = file.path(out_dir, paste0("module_", tag, ".graphml")),
                        format = "graphml")
  }
  
  log_msg("---- ", out_prefix, " | PPI=", ppi_name, " | done ----", logfile = logfile)
  invisible(TRUE)
}

# Run across your 3 PPIs: INTEGRATED proteomics + GEX
fdr_thresholds <- c(0.05)

# recommended: TRUE to avoid BUM fit dominated by p=1 nodes
restrict_flag <- TRUE

total_steps <- length(names(ppi_list)) * length(fdr_thresholds)
step <- 0L

for (ppi_name in names(ppi_list)) {
  log_msg("=== PPI: ", ppi_name, " (integrated run) ===", logfile = logfile_global)
  
  step <- step + length(fdr_thresholds)
  log_msg("Progress: starting integrated run for ", ppi_name,
          " | approx step ", step - length(fdr_thresholds) + 1, "-", step, "/", total_steps,
          logfile = logfile_global)
  
  run_bionet_integrated_ppi(ppi_list[[ppi_name]], ppi_name, omics,
                            fdr_vec = fdr_thresholds,
                            keep_lcc = TRUE,
                            restrict_to_measured = restrict_flag,
                            measured_mode = "any",  # change to "both" if you require both omics measured
                            logfile = logfile_global)
}

log_msg("All done. Results in output/bionet_integrated_*; log: output/bionet_runlog_integrated.txt",
        logfile = logfile_global)