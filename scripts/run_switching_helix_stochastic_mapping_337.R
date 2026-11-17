#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ape)
  library(phytools)
})

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(2, 3)) {
  stop("Usage: run_switching_helix_stochastic_mapping_337.R <nsim> <run_label> [ER|SYM|ARD]")
}
nsim <- as.integer(args[[1]])
run_label <- args[[2]]
model_override <- if (length(args) == 3) args[[3]] else NA_character_
if (is.na(nsim) || nsim < 1) stop("nsim must be a positive integer")

project <- normalizePath("Vertebrate_Gene_History_Redesign_20261102")
root <- file.path(project, "analysis/switching_helix_asr_337")
tree_path <- file.path(root, "runs/TAXONOMY_STRICT_LG_F_R4/DNMT3L_289_taxonomy_strict.treefile")
states_path <- file.path(project, "results/switching_helix_337_observed_states.tsv")
out_dir <- file.path(root, "runs", paste0("STOCHASTIC_MAP_NP349_", run_label))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

tree <- read.tree(tree_path)
states_table <- read.delim(states_path, stringsAsFactors = FALSE, check.names = FALSE)
site <- subset(states_table, human_np_position == 349)
site$strict_cross_method_concordance <- as.logical(site$strict_cross_method_concordance)
site$state3 <- ifelse(
  !site$strict_cross_method_concordance | site$mafft_human_coordinate_state %in% c("-", "X", "?"),
  NA,
  ifelse(site$mafft_human_coordinate_state == "Q", "Q",
    ifelse(site$mafft_human_coordinate_state == "K", "K", "O")
  )
)
x <- setNames(site$state3, site$safe_id)
x <- x[tree$tip.label]
if (any(is.na(names(x))) || !identical(names(x), tree$tip.label)) stop("Tip-state ordering failed")
missing_tips <- names(x)[is.na(x)]
if (length(missing_tips)) {
  tree <- drop.tip(tree, missing_tips)
  x <- x[tree$tip.label]
}
x <- factor(x, levels = c("Q", "K", "O"))
names(x) <- tree$tip.label

amphibian_tips <- intersect(site$safe_id[site$clade == "Amphibia" & !is.na(site$state3)], tree$tip.label)
if (length(amphibian_tips) < 2) stop("Insufficient amphibian outgroup tips")
if (!is.monophyletic(tree, amphibian_tips)) stop("Amphibian outgroup is not monophyletic")
tree <- root(tree, outgroup = amphibian_tips, resolve.root = TRUE)
if (!is.rooted(tree)) stop("Tree rooting failed")

models <- c("ER", "SYM", "ARD")
fits <- lapply(models, function(model) fitMk(tree, x, model = model, pi = "fitzjohn"))
names(fits) <- models
n <- length(x)
model_rows <- lapply(models, function(model) {
  fit <- fits[[model]]
  k <- length(fit$rates)
  aic <- 2 * k - 2 * fit$logLik
  aicc <- aic + (2 * k * (k + 1)) / (n - k - 1)
  data.frame(model = model, tips = n, states = 3, parameters = k,
             log_likelihood = fit$logLik, AIC = aic, AICc = aicc)
})
model_table <- do.call(rbind, model_rows)
model_table$delta_AICc <- model_table$AICc - min(model_table$AICc)
model_table$akaike_weight <- exp(-0.5 * model_table$delta_AICc) /
  sum(exp(-0.5 * model_table$delta_AICc))
model_table <- model_table[order(model_table$AICc), ]
best_model <- model_table$model[[1]]
if (!is.na(model_override)) {
  if (!model_override %in% models) stop("Model override must be ER, SYM, or ARD")
  best_model <- model_override
}
write.table(model_table, file.path(out_dir, "NP349_transition_model_comparison.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.tree(tree, file.path(out_dir, "DNMT3L_282_resolved_NP349_rooted_taxonomy_tree.nwk"))
write.table(data.frame(tip = missing_tips), file.path(out_dir, "NP349_masked_tips.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

set.seed(20261102)
maps <- make.simmap(tree, x, model = best_model, nsim = nsim,
                    pi = "fitzjohn", message = FALSE)
if (nsim == 1) maps <- list(maps)
class(maps) <- c("multiSimmap", "multiPhylo")
saveRDS(maps, file.path(out_dir, "NP349_stochastic_maps.rds"), compress = "xz")

state_order <- c("Q", "K", "O")
map_rows <- vector("list", length(maps))
edge_accumulator <- vector("list", nrow(tree$edge))
for (edge_i in seq_len(nrow(tree$edge))) {
  edge_accumulator[[edge_i]] <- c(any = 0, into_Q = 0, into_K = 0, direct_QK = 0)
}

node_label <- function(node) {
  if (node <= length(tree$tip.label)) return(tree$tip.label[[node]])
  label <- tree$node.label[[node - length(tree$tip.label)]]
  if (is.null(label) || is.na(label) || label == "") paste0("Node", node) else label
}

for (i in seq_along(maps)) {
  counted <- countSimmap(maps[[i]], states = state_order, message = FALSE)$Tr
  into_q <- sum(counted[c("K", "O"), "Q"])
  into_k <- sum(counted[c("Q", "O"), "K"])
  map_rows[[i]] <- data.frame(
    map = i,
    total_transitions = sum(counted),
    transitions_into_Q = into_q,
    transitions_into_K = into_k,
    direct_Q_to_K = counted["Q", "K"],
    direct_K_to_Q = counted["K", "Q"],
    recurrent_Q_gain = into_q >= 2,
    recurrent_K_gain = into_k >= 2
  )
  for (edge_i in seq_len(nrow(tree$edge))) {
    history <- names(maps[[i]]$maps[[edge_i]])
    if (length(history) < 2) next
    from <- history[-length(history)]
    to <- history[-1]
    edge_accumulator[[edge_i]]["any"] <- edge_accumulator[[edge_i]]["any"] + 1
    edge_accumulator[[edge_i]]["into_Q"] <- edge_accumulator[[edge_i]]["into_Q"] + any(to == "Q" & from != "Q")
    edge_accumulator[[edge_i]]["into_K"] <- edge_accumulator[[edge_i]]["into_K"] + any(to == "K" & from != "K")
    edge_accumulator[[edge_i]]["direct_QK"] <- edge_accumulator[[edge_i]]["direct_QK"] +
      any((from == "Q" & to == "K") | (from == "K" & to == "Q"))
  }
}
map_table <- do.call(rbind, map_rows)
write.table(map_table, file.path(out_dir, "NP349_transition_counts_by_map.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

quant <- function(values, probability) unname(quantile(values, probability, type = 8))
metrics <- c("total_transitions", "transitions_into_Q", "transitions_into_K",
             "direct_Q_to_K", "direct_K_to_Q")
summary_rows <- lapply(metrics, function(metric) {
  values <- map_table[[metric]]
  data.frame(metric = metric, nsim = nsim, mean = mean(values), median = median(values),
             q025 = quant(values, 0.025), q975 = quant(values, 0.975),
             minimum = min(values), maximum = max(values))
})
summary_table <- do.call(rbind, summary_rows)
summary_table <- rbind(summary_table,
  data.frame(metric = "probability_recurrent_Q_gain", nsim = nsim,
             mean = mean(map_table$recurrent_Q_gain), median = NA, q025 = NA, q975 = NA,
             minimum = NA, maximum = NA),
  data.frame(metric = "probability_recurrent_K_gain", nsim = nsim,
             mean = mean(map_table$recurrent_K_gain), median = NA, q025 = NA, q975 = NA,
             minimum = NA, maximum = NA)
)
write.table(summary_table, file.path(out_dir, "NP349_stochastic_mapping_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

edge_rows <- lapply(seq_len(nrow(tree$edge)), function(edge_i) {
  parent <- tree$edge[edge_i, 1]
  child <- tree$edge[edge_i, 2]
  values <- edge_accumulator[[edge_i]] / nsim
  child_tip_count <- if (child <= length(tree$tip.label)) 1 else length(extract.clade(tree, child)$tip.label)
  data.frame(edge = edge_i, parent_node = node_label(parent), child_node = node_label(child),
             child_tip_count = child_tip_count,
             probability_any_transition = values[["any"]],
             probability_transition_into_Q = values[["into_Q"]],
             probability_transition_into_K = values[["into_K"]],
             probability_direct_QK_transition = values[["direct_QK"]])
})
edge_table <- do.call(rbind, edge_rows)
edge_table <- edge_table[order(-edge_table$probability_any_transition), ]
write.table(edge_table, file.path(out_dir, "NP349_branch_transition_probabilities.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("tips=", n, " masked=", length(missing_tips), " mapping_model=", best_model,
    " nsim=", nsim, "\n", sep = "")
print(model_table)
print(summary_table)
