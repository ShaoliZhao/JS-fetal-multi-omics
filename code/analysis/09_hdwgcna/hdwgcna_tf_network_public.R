# TF-focused hdWGCNA network summaries.

source("analysis/00_setup/packages.R")
source("analysis/00_setup/project_config.R")

suppressPackageStartupMessages({
  library(hdWGCNA)
  library(igraph)
  library(ggraph)
})

out_root <- file.path(paths$results, "09_hdwgcna_tf_networks")
fig_root <- file.path(paths$figures, "09_hdwgcna_tf_networks")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_root, recursive = TRUE, showWarnings = FALSE)

datasets <- list(
  cerebellum = file.path(paths$results, "09_hdwgcna", "cerebellum", "cerebellum_hdWGCNA.rds"),
  kidney = file.path(paths$results, "09_hdwgcna", "kidney", "kidney_hdWGCNA.rds")
)

tf_file <- file.path(paths$data_external, "gene_sets", "human_transcription_factors.txt")
if (file.exists(tf_file)) {
  tf_genes <- unique(readLines(tf_file))
} else {
  tf_genes <- c("ATOH1", "PAX6", "EOMES", "SOX2", "HES1", "HES5", "FOXC1", "FOXD1", "PBX1", "TCF7L2", "YAP1")
}

for (dataset in names(datasets)) {
  if (!file.exists(datasets[[dataset]])) {
    message("Skip ", dataset, ": hdWGCNA object not found at ", datasets[[dataset]])
    next
  }

  out_dir <- file.path(out_root, dataset)
  fig_dir <- file.path(fig_root, dataset)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

  obj <- readRDS(datasets[[dataset]])
  modules <- GetModules(obj) %>% filter(module != "grey")
  hub_df <- GetHubGenes(obj, n_hubs = 50)

  tf_modules <- modules %>% filter(gene_name %in% tf_genes)
  tf_hubs <- hub_df %>% filter(gene_name %in% tf_genes)
  write.csv(tf_modules, file.path(out_dir, "module_transcription_factors.csv"), row.names = FALSE)
  write.csv(tf_hubs, file.path(out_dir, "hub_transcription_factors.csv"), row.names = FALSE)

  edges <- modules %>%
    filter(gene_name %in% tf_genes | kME > quantile(kME, 0.95, na.rm = TRUE)) %>%
    group_by(module) %>%
    arrange(desc(kME), .by_group = TRUE) %>%
    slice_head(n = 25) %>%
    ungroup() %>%
    select(module, gene_name, kME)

  write.csv(edges, file.path(out_dir, "top_module_network_nodes.csv"), row.names = FALSE)

  graph_edges <- edges %>%
    mutate(from = module, to = gene_name, weight = kME) %>%
    select(from, to, weight)
  g <- graph_from_data_frame(graph_edges, directed = FALSE)

  p <- ggraph(g, layout = "fr") +
    geom_edge_link(aes(alpha = weight), color = "grey70") +
    geom_node_point(aes(color = grepl("^M", name), size = ifelse(name %in% tf_genes, 3, 1.5))) +
    geom_node_text(aes(label = ifelse(name %in% tf_genes | grepl("^M", name), name, "")), repel = TRUE, size = 3) +
    scale_color_manual(values = c("TRUE" = "#3A6EA5", "FALSE" = "#B23A48")) +
    theme_void() +
    theme(legend.position = "none")
  ggsave(file.path(fig_dir, "tf_module_network.pdf"), p, width = 8, height = 6)
}
