organism_to_orgdb <- function(organism) {
  switch(
    tolower(organism),
    human = "org.Hs.eg.db",
    mouse = "org.Mm.eg.db",
    rat = "org.Rn.eg.db",
    "org.Hs.eg.db"
  )
}

organism_to_kegg <- function(organism) {
  switch(
    tolower(organism),
    human = "hsa",
    mouse = "mmu",
    rat = "rno",
    "hsa"
  )
}

resolve_universe <- function(config, background_genes, gene_info, de_tbl = NULL) {
  if (config$enrichment$universe == "custom" && !is.null(background_genes)) {
    return(background_genes)
  }
  if (config$enrichment$universe == "expressed_genes" && !is.null(de_tbl)) {
    return(de_tbl$gene_id)
  }
  gene_info$gene_id
}

safe_bitr <- function(ids, from_type, to_type, orgdb) {
  tryCatch(
    clusterProfiler::bitr(ids, fromType = from_type, toType = to_type, OrgDb = orgdb),
    error = function(e) NULL
  )
}

run_ora_for_gene_set <- function(genes, config, universe, label, out_dir) {
  if (length(genes) < 5 || !isTRUE(config$enrichment$enabled)) return(NULL)
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) return(NULL)
  orgdb <- organism_to_orgdb(config$input$organism)
  key_type <- config$input$id_type %||% "SYMBOL"
  mapped <- safe_bitr(genes, key_type, "ENTREZID", orgdb)
  if (is.null(mapped) || nrow(mapped) == 0) return(NULL)
  gene_ids <- unique(mapped$ENTREZID)
  universe_map <- safe_bitr(universe, key_type, "ENTREZID", orgdb)
  universe_ids <- if (is.null(universe_map)) NULL else unique(universe_map$ENTREZID)

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  outputs <- list()

  if ("GO_BP" %in% config$enrichment$databases) {
    ego <- tryCatch(
      clusterProfiler::enrichGO(
        gene = gene_ids,
        OrgDb = orgdb,
        keyType = "ENTREZID",
        ont = "BP",
        universe = universe_ids,
        readable = TRUE
      ),
      error = function(e) NULL
    )
    if (!is.null(ego)) {
      readr::write_csv(as.data.frame(ego), file.path(out_dir, paste0(label, "_go_bp.csv")))
      outputs$go_bp <- file.path(out_dir, paste0(label, "_go_bp.csv"))
    }
  }

  if ("KEGG" %in% config$enrichment$databases) {
    ekegg <- tryCatch(
      clusterProfiler::enrichKEGG(
        gene = gene_ids,
        organism = organism_to_kegg(config$input$organism),
        keyType = "ncbi-geneid",
        universe = universe_ids
      ),
      error = function(e) NULL
    )
    if (!is.null(ekegg)) {
      readr::write_csv(as.data.frame(ekegg), file.path(out_dir, paste0(label, "_kegg.csv")))
      outputs$kegg <- file.path(out_dir, paste0(label, "_kegg.csv"))
    }
  }

  outputs
}
