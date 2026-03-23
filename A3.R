library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(grid)
library(vegan)
library(ALDEx3)

# Metadata 
metadata <- read.csv("metadata.csv",stringsAsFactors = FALSE,check.names = FALSE)

# Combined species abundance matrix (samples = columns, taxa = rows)
bracken_combined <- read.delim("bracken_combined_abundance.tsv",row.names = 1,check.names = FALSE,stringsAsFactors = FALSE)

# numeric matrix for downstream stats (phyloseq, vegan, etc.)
bracken_combined_mat <- as.matrix(bracken_combined)
mode(bracken_combined_mat) <- "numeric"

bracken_files <- list.files(pattern = "^bracken_SRR[0-9]+\\.txt$",full.names = TRUE)

bracken_per_sample <- lapply(bracken_files, function(f) {read.delim(f, stringsAsFactors = FALSE, check.names = FALSE)})
names(bracken_per_sample) <- sub("^bracken_(SRR[0-9]+)\\.txt$","\\1",basename(bracken_files))

# Quick sanity check: metadata samples should match abundance columns
stopifnot(setequal(metadata$sample, colnames(bracken_combined)))

# taxonomic abundance barplot
# Long format: one row per sample × taxon
abund_long <- bracken_combined_mat %>%
  as.data.frame(check.names = FALSE) %>%
  rownames_to_column("name") %>%
  pivot_longer(-name, names_to = "sample", values_to = "abundance")

meta <- metadata %>%
  mutate(group = factor(group, levels = c("omnivore", "vegan")),sample = as.character(sample))

plot_df <- abund_long %>%
  inner_join(meta, by = "sample")

# Top 20 taxa by total abundance across all samples
top20 <- plot_df %>%
  group_by(name) %>%
  summarise(total = sum(abundance), .groups = "drop") %>%
  arrange(desc(total)) %>%
  slice(1:20) %>%
  pull(name)

# Legend order: most abundant taxon first (matches reference figure)
plot_df <- plot_df %>%
  filter(name %in% top20) %>%
  mutate(name = factor(name, levels = top20))

# X order: omnivore samples then vegan (left-to-right in reference)
sample_levels <- c(
  "SRR8146935", "SRR8146936", "SRR8146938",
  "SRR8146963", "SRR8146968", "SRR8146973")
plot_df <- plot_df %>%
  mutate(sample = factor(sample, levels = sample_levels))

stack_max <- plot_df %>%
  group_by(sample) %>%
  summarise(s = sum(abundance), .groups = "drop") %>%
  pull(s) %>%
  max()
y_upper <- max(0.8, stack_max * 1.02)
y_breaks <- seq(0, ceiling(y_upper / 0.2) * 0.2, by = 0.2)

# Discrete fill palette (broad spectrum, similar to reference)
n_col <- length(levels(plot_df$name))
fill_cols <- hcl.colors(n_col, palette = "Spectral")
names(fill_cols) <- levels(plot_df$name)

p <- ggplot(plot_df, aes(x = sample, y = abundance, fill = name)) +
  geom_col(width = 0.9, colour = NA) +
  facet_wrap(~group, nrow = 1, scales = "free_x") +
  scale_y_continuous(
    name = "Relative Abundance",
    limits = c(0, y_upper),
    breaks = y_breaks,
    expand = c(0, 0) ) +
  scale_x_discrete(name = "Sample") +
  scale_fill_manual(
    name = "name",
    values = fill_cols,
    breaks = top20) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey85", colour = "black"),
    strip.text = element_text(colour = "black"),
    legend.text = element_text(face = "italic"),
    legend.key.size = unit(0.45, "cm"),
    legend.title = element_text(face = "plain"),
    panel.border = element_rect(colour = "black", fill = NA))

ggsave(
  "taxonomic_abundance_by_group.png",
  plot = p,
  width = 10,
  height = 5,
  dpi = 300,
  bg = "white")

message("Saved: taxonomic_abundance_by_group.png")

#alpha diversity metrics
# build integer count matrix: rows = taxa, cols = samples (required by vegan)
samples <- as.character(metadata$sample)
stopifnot(setequal(names(bracken_per_sample), samples))

taxa <- sort(unique(unlist(lapply(bracken_per_sample, function(df) df$name))))
count_mat <- matrix(
  0L,
  nrow = length(taxa),
  ncol = length(samples),
  dimnames = list(taxa, samples))
for (s in samples) {df <- bracken_per_sample[[s]]
  idx <- match(df$name, taxa)
  ok <- !is.na(idx)
  count_mat[idx[ok], s] <- as.integer(round(df$new_est_reads[ok]))}

ns <- length(samples)
if (ncol(count_mat) != ns && nrow(count_mat) == ns) {count_mat <- t(count_mat)}
if (ncol(count_mat) != ns) {stop("count_mat must have one column per sample (ncol = ", ns, "). ",
    "Got ", nrow(count_mat), " x ", ncol(count_mat), ".")}
colnames(count_mat) <- samples
sample_names <- colnames(count_mat)

# Chao1 + observed richness (vegan::estimateR)
est_r <- sapply(seq_len(ncol(count_mat)), function(i) {
  vegan::estimateR(count_mat[, i])})
if (is.null(rownames(est_r))) {
  rownames(est_r) <- names(vegan::estimateR(count_mat[, 1]))}
colnames(est_r) <- sample_names

chao_row <- grep("^S\\.chao1$", rownames(est_r), value = TRUE)[1]
if (is.na(chao_row)) {
  chao_row <- grep("chao1", rownames(est_r), ignore.case = TRUE, value = TRUE)[1]}
obs_row <- grep("^S\\.obs$", rownames(est_r), value = TRUE)[1]
if (is.na(obs_row)) {
  obs_row <- grep("^S\\.obs", rownames(est_r), value = TRUE)[1]}
if (is.na(obs_row)) {obs_row <- rownames(est_r)[1]}

if (is.na(chao_row)) {
  stop("Could not find Chao1 row in estimateR() output; check vegan version.")}

chao1_df <- tibble(
  sample = colnames(est_r),
  chao1 = as.numeric(est_r[chao_row, , drop = TRUE]),
  s_obs = as.numeric(est_r[obs_row, , drop = TRUE])) %>%
  left_join(
    metadata %>% mutate(sample = as.character(sample)),
    by = "sample") %>%
  mutate(
    group = factor(group, levels = c("omnivore", "vegan")),
    sample = factor(
      sample,
      levels = c(
        "SRR8146935", "SRR8146936", "SRR8146938",
        "SRR8146963", "SRR8146968", "SRR8146973")))

if (anyNA(chao1_df$sample) || anyNA(chao1_df$group)) {
  stop(
    "Sample/group lost after join. colnames(est_r): ",
    paste(colnames(est_r), collapse = ", "))}

write.csv(chao1_df, "chao1_per_sample.csv", row.names = FALSE)

# visualization: Chao1 per sample, filled by group
p <- ggplot(chao1_df, aes(x = sample, y = chao1, fill = group)) +
  geom_col(width = 0.75, colour = "grey30") +
  scale_fill_manual(
    values = c(omnivore = "#6BAED6", vegan = "#74C476"),
    name = "Group") +
  labs(
    x = "Sample",
    y = "Chao1 estimated richness",
    title = "Chao1 diversity (species-level Bracken counts)") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid.minor = element_blank())

ggsave(
  "chao1_barplot.png",
  plot = p,
  width = 8,
  height = 5,
  dpi = 300,
  bg = "white")

message("Wrote: chao1_per_sample.csv, chao1_barplot.png")

#berger-parker dominance/evenness metric

samples <- as.character(metadata$sample)
stopifnot(setequal(colnames(bracken_combined_mat), samples))

# berger-parker dominance: max abundance / total (same for counts or proportions)
berger_parker <- apply(bracken_combined_mat, 2, function(x) {
  x <- as.numeric(x)
  sum_x <- sum(x)
  if (sum_x <= 0) {
    return(NA_real_)}
  max(x) / sum_x})

bp_df <- tibble(
  sample = names(berger_parker),
  berger_parker_dominance = as.numeric(berger_parker)) %>%
  left_join(metadata %>% mutate(sample = as.character(sample)), by = "sample") %>%
  mutate(group = factor(group, levels = c("omnivore", "vegan")),
    sample = factor(
      sample,
      levels = c("SRR8146935", "SRR8146936", "SRR8146938",
        "SRR8146963", "SRR8146968", "SRR8146973")))

write.csv(bp_df, "berger_parker_per_sample.csv", row.names = FALSE)

p <- ggplot(bp_df, aes(x = sample, y = berger_parker_dominance, fill = group)) +
  geom_col(width = 0.75, colour = "grey30") +
  scale_fill_manual(
    values = c(omnivore = "#6BAED6", vegan = "#74C476"),
    name = "Group") +
  scale_y_continuous(limits = c(0, NA), expand = c(0, 0)) +
  labs(x = "Sample",
    y = "Berger-Parker dominance (D)",
    title = "Berger-Parker dominance by sample",
    subtitle = "Share of reads in the most abundant species; lower = more even") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid.minor = element_blank())

ggsave(
  "berger_parker_barplot.png",
  plot = p,
  width = 8,
  height = 5,
  dpi = 300,
  bg = "white")

message("Wrote: berger_parker_per_sample.csv, berger_parker_barplot.png")

#shannon diversity index

samples <- as.character(metadata$sample)
mat <- bracken_combined_mat
mode(mat) <- "numeric"

# rows = taxa, columns = samples
ns <- length(samples)
if (ncol(mat) != ns && nrow(mat) == ns) {
  mat <- t(mat)}
if (ncol(mat) != ns) {
  stop(
    "Abundance matrix must have one column per sample (ncol = ", ns, "). ",
    "Got ", nrow(mat), " x ", ncol(mat), ".")}
if (!all(samples %in% colnames(mat))) {
  stop(
    "Sample IDs from metadata are not column names of the abundance matrix.")}
mat <- mat[, samples, drop = FALSE]
colnames(mat) <- samples

# One H' per sample column
shannon <- sapply(seq_len(ncol(mat)), function(i) {
  vegan::diversity(mat[, i], index = "shannon")})
names(shannon) <- colnames(mat)
shannon_df <- tibble(
  sample = names(shannon),
  shannon = as.numeric(shannon)) %>%
  left_join(metadata %>% mutate(sample = as.character(sample)), by = "sample") %>%
  mutate(
    group = factor(group, levels = c("omnivore", "vegan")),
    sample = factor(
      sample,
      levels = c(
        "SRR8146935", "SRR8146936", "SRR8146938",
        "SRR8146963", "SRR8146968", "SRR8146973")))

if (anyNA(shannon_df$sample) || anyNA(shannon_df$group)) {
  stop("Sample or group is NA after join; check metadata vs colnames(mat).")}
if (nrow(shannon_df) != ns) {
  stop("Expected ", ns, " Shannon values, got ", nrow(shannon_df), ".")}

write.csv(shannon_df, "shannon_per_sample.csv", row.names = FALSE)

p <- ggplot(shannon_df, aes(x = sample, y = shannon, fill = group)) +
  geom_col(width = 0.75, colour = "grey30") +
  scale_fill_manual(
    values = c(omnivore = "#6BAED6", vegan = "#74C476"),
    name = "Group") +
  scale_y_continuous(limits = c(0, NA), expand = c(0, 0)) +
  labs(
    x = "Sample",
    y = "Shannon diversity H' (ln; nats)",
    title = "Shannon diversity by sample",
    subtitle = "Natural log (ln); Bracken species relative abundances") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid.minor = element_blank())

ggsave(
  "shannon_barplot.png",
  plot = p,
  width = 8,
  height = 5,
  dpi = 300,
  bg = "white")

message("Wrote: shannon_per_sample.csv, shannon_barplot.png")

#bray-curtis beta diversity metric
samples <- as.character(metadata$sample)
mat <- bracken_combined_mat
mode(mat) <- "numeric"

ns <- length(samples)
if (ncol(mat) != ns && nrow(mat) == ns) {
  mat <- t(mat)}
if (ncol(mat) != ns) {
  stop(
    "Abundance matrix must have one column per sample (ncol = ", ns, "). ",
    "Got ", nrow(mat), " x ", ncol(mat), ".")}
if (!all(samples %in% colnames(mat))) {
  stop("Sample IDs from metadata are not column names of the abundance matrix.")}
mat <- mat[, samples, drop = FALSE]
mat <- as.matrix(mat)
colnames(mat) <- samples
stopifnot(ncol(mat) == ns, length(colnames(mat)) == ns)

mat[is.na(mat)] <- 0
mat[is.infinite(mat)] <- 0
mat[mat < 0] <- 0

cs <- colSums(mat)
if (any(cs <= 0 | !is.finite(cs))) {
  stop(
    "Samples with zero or invalid total abundance: ",
    paste(colnames(mat)[cs <= 0 | !is.finite(cs)], collapse = ", "))}

#bray-curtis (between samples = rows of comm)
comm <- t(mat)
rownames(comm) <- colnames(mat)
colnames(comm) <- rownames(mat)

bray_dist <- vegdist(comm, method = "bray")
bray_mat <- as.matrix(bray_dist)

if (anyNA(bray_mat) || !all(is.finite(bray_mat))) {
  # Tiny pseudocount avoids rare 0/0 cases in Bray–Curtis; negligible for proportions
  mat <- mat + 1e-10
  comm <- t(mat)
  rownames(comm) <- colnames(mat)
  colnames(comm) <- rownames(mat)
  bray_dist <- vegdist(comm, method = "bray")
  bray_mat <- as.matrix(bray_dist)}
if (anyNA(bray_mat)) {
  stop("Bray-Curtis distances still contain NA after cleaning abundances. ",
    "Inspect bracken_combined_mat for bad values.")}

if (nrow(bray_mat) != ns || ncol(bray_mat) != ns) {
  stop(
    "Bray-Curtis matrix is ", nrow(bray_mat), " x ", ncol(bray_mat),
    " but expected ", ns, " x ", ns,
    " (ncol(mat) = ", ncol(mat), ", nrow(mat) = ", nrow(mat), ").")}
dimnames(bray_mat) <- list(samples, samples)

write.csv(cbind(sample = rownames(bray_mat), as.data.frame(bray_mat)),
  "bray_curtis_matrix.csv",
  row.names = FALSE)

# PCoA (classical scaling) on Bray-Curtis
pco <- cmdscale(bray_mat, k = min(5, ns - 1), eig = TRUE)
rownames(pco$points) <- samples

pco_df <- tibble(
  sample = rownames(pco$points),
  Axis1 = pco$points[, 1],
  Axis2 = pco$points[, 2]) %>%
  mutate(sample = as.character(sample)) %>%
  left_join(metadata %>% mutate(sample = as.character(sample)), by = "sample") %>%
  mutate(
    group = factor(group, levels = c("omnivore", "vegan")),
    sample = factor(sample, levels = samples))

if (anyNA(pco_df$group)) {
  stop(
    "Metadata join failed for PCoA (group is NA). Samples: ",
    paste(pco_df$sample, collapse = ", "))}

eig_pos <- pco$eig[Re(pco$eig) > 0]
var_exp <- eig_pos / sum(eig_pos)
lab1 <- paste0("PCoA1 (", format(100 * var_exp[1], digits = 1), "%)")
lab2 <- paste0("PCoA2 (", format(100 * var_exp[2], digits = 1), "%)")

p_pcoa <- ggplot(pco_df, aes(x = Axis1, y = Axis2, colour = group, label = sample)) +
  geom_point(size = 4, alpha = 0.9) +
  geom_text(vjust = -0.8, size = 2.8, show.legend = FALSE, colour = "grey25") +
  scale_colour_manual(
    values = c(omnivore = "#6BAED6", vegan = "#74C476"),
    name = "Group",
    drop = FALSE) +
  labs(
    x = lab1,
    y = lab2,
    title = "PCoA on Bray-Curtis dissimilarity",
    subtitle = "Bracken species relative abundances") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())

ggsave(
  "bray_curtis_pcoa.png",
  plot = p_pcoa,
  width = 7,
  height = 5,
  dpi = 300,
  bg = "white")

# summary statistics (console + bray_curtis_stats.txt)
grp <- as.character(metadata$group[match(samples, metadata$sample)])
ut <- upper.tri(bray_mat)
same_grp <- outer(grp, grp, function(a, b) a == b)
off <- bray_mat[ut]
within_d <- bray_mat[ut & same_grp[ut]]
between_d <- bray_mat[ut & !same_grp[ut]]

stats_lines <- capture.output({
  cat("========== Bray-Curtis / PCoA summary ==========\n\n")
  cat("Full dissimilarity matrix (rounded):\n")
  print(round(bray_mat, 4))
  cat("\nAll pairwise distances (upper triangle, excluding diagonal):\n")
  print(summary(off))
  cat(sprintf("\nMean pairwise distance: %.4f\n", mean(off)))
  if (length(within_d)) {
    cat(sprintf("Mean within-group (same diet): %.4f\n", mean(within_d)))} else {
    cat("Mean within-group: (no same-diet pairs in upper triangle)\n")}
  if (length(between_d)) {
    cat(sprintf("Mean between-group (different diet): %.4f\n", mean(between_d)))} else {
    cat("Mean between-group: (no cross-diet pairs)\n")}
  cat("\n--- PCoA (proportion of variance, positive eigenvalues) ---\n")
  for (i in seq_len(min(5L, length(var_exp)))) {
    cat(sprintf("  PCoA axis %i: %.2f%%\n", i, 100 * Re(var_exp[i])))}
  cat("\nPCoA site scores (first two axes):\n")
  print(select(pco_df, sample, Axis1, Axis2) %>%
      mutate(Axis1 = round(Axis1, 4), Axis2 = round(Axis2, 4)))})

writeLines(stats_lines, "bray_curtis_stats.txt")
cat(paste(stats_lines, collapse = "\n"), "\n")

message("Wrote: bray_curtis_matrix.csv, bray_curtis_heatmap.png, bray_curtis_pcoa.png, ",
  "bray_curtis_stats.txt (see console output for the same stats)")

if (!requireNamespace("ALDEx3", quietly = TRUE)) {
  stop("Install ALDEx3: install.packages(\"ALDEx3\")", call. = FALSE)}

samples <- as.character(metadata$sample)
stopifnot(setequal(names(bracken_per_sample), samples))

# integer count matrix Y: rows = taxa (D), columns = samples (N) — required by ALDEx3
taxa <- sort(unique(unlist(lapply(bracken_per_sample, function(df) df$name))))
Y <- matrix(0L,
  nrow = length(taxa),
  ncol = length(samples),
  dimnames = list(taxa, samples))
for (s in samples) {
  df <- bracken_per_sample[[s]]
  idx <- match(df$name, taxa)
  ok <- !is.na(idx)
  Y[idx[ok], s] <- as.integer(round(df$new_est_reads[ok]))}

Y <- Y[, samples, drop = FALSE]

# drop taxa with almost no reads
min_total <- 10L
min_samples <- 2L
keep <- rowSums(Y) >= min_total & rowSums(Y > 0L) >= min_samples
Y <- Y[keep, , drop = FALSE]
if (nrow(Y) < 1L) {
  stop("No taxa left after filtering; lower min_total or min_samples.")}

# Covariates: rows must match column order of Y
design <- data.frame(
  group = factor(
    metadata$group[match(colnames(Y), metadata$sample)],
    levels = c("omnivore", "vegan")),
  row.names = colnames(Y))

# ALDEx3
set.seed(1)
res <- aldex(
  Y,
  ~ group,
  design,
  method = "lm",
  nsample = 1000L,
  scale = clr.sm,
  gamma = 0.5,
  test = "t.HC3",
  p.adjust.method = "BH",
  n.cores = 1L)

rn <- rownames(res$p.val)
if (is.null(rn)) {
  coef_row <- 2L} else {
  coef_row <- which(rn == "groupvegan")[1]
  if (is.na(coef_row)) {
    coef_row <- grep("vegan", rn, ignore.case = TRUE)[1]}
  if (is.na(coef_row)) {
    coef_row <- nrow(res$p.val)}}

est_mean <- apply(res$estimate, c(1L, 2L), mean)
se_mean <- apply(res$std.error, c(1L, 2L), mean)

out <- data.frame(
  taxon = rownames(Y),
  log_ratio_effect = as.numeric(est_mean[coef_row, , drop = TRUE]),
  se_mean = as.numeric(se_mean[coef_row, , drop = TRUE]),
  p_value = as.numeric(res$p.val[coef_row, , drop = TRUE]),
  p_adj_BH = as.numeric(res$p.val.adj[coef_row, , drop = TRUE]),
  row.names = NULL,
  stringsAsFactors = FALSE)
out <- out[order(out$p_adj_BH, out$p_value), , drop = FALSE]

write.csv(out, "aldex3_vegan_vs_omnivore.csv", row.names = FALSE)

#visualizations (ALDEx3 log-ratio scale; omnivore = reference) ---
plot_df <- out %>%
  mutate(p_adj_plot = ifelse(is.na(p_adj_BH), 1, p_adj_BH),
    neglog10_padj = -log10(pmax(p_adj_plot, 1e-300)),
    significant = if_else(!is.na(p_adj_BH) & p_adj_BH < 0.05, "FDR < 0.05", "FDR >= 0.05"))

n_bar <- min(25L, nrow(out))
bar_df <- slice_head(out, n = n_bar)

p_bar <- ggplot(bar_df, aes(x = reorder(taxon, log_ratio_effect), y = log_ratio_effect)) +
  geom_col(aes(fill = log_ratio_effect > 0), width = 0.85) +
  geom_errorbar(aes(
      ymin = log_ratio_effect - 1.96 * se_mean,
      ymax = log_ratio_effect + 1.96 * se_mean),
    width = 0.2,
    linewidth = 0.25) +
  scale_fill_manual(
    values = c(`FALSE` = "#6BAED6", `TRUE` = "#74C476"),
    breaks = c(FALSE, TRUE),
    labels = c("Higher in omnivore", "Higher in vegan"),
    name = NULL) +
  coord_flip() +
  labs(
    title = paste0("Top ", n_bar, " taxa by smallest adjusted p-value"),
    subtitle = "Bars = mean ALDEx3 coefficient; error bars ~95% CI (1.96 * SE)",
    x = NULL,
    y = "Effect (vegan vs omnivore)") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 6))

ggsave(
  "aldex3_top_taxa.png",
  plot = p_bar,
  width = 9,
  height = min(14, max(5, n_bar * 0.22)),
  dpi = 300,
  bg = "white")

cat("\n=== ALDEx3: vegan vs omnivore (reference: omnivore) ===\n")
cat("Taxa analyzed:", nrow(Y), "| Samples:", ncol(Y), "\n")
cat("Coefficient row used:", coef_row, if (!is.null(rn)) paste0(" (", rn[coef_row], ")") else "", "\n")
cat("Top 15 by adjusted p-value:\n")
print(head(out, 15L))
cat("\nFull table: aldex3_vegan_vs_omnivore.csv\n")
cat("Interpretation: positive log-ratio effect ≈ higher relative abundance in vegan vs omnivore.\n")

message("Wrote: aldex3_vegan_vs_omnivore.csv, aldex3_top_taxa.png")
