# Shared pieces for the ORA and GSEA modules -------------------------------------
#
# Pathway names are borrowed from well-known collections for realism; their gene
# memberships here are simulated.

pathway_names <- c(
  "Interferon alpha response", "Interferon gamma response", "Inflammatory response", "TNF-alpha signaling via NF-kB",
  "E2F targets (cell cycle)", "G2/M checkpoint", "Mitotic spindle", "DNA repair",
  "Oxidative phosphorylation", "Glycolysis", "Fatty acid metabolism", "Cholesterol homeostasis",
  "Hypoxia", "Apoptosis", "p53 pathway", "Unfolded protein response",
  "MYC targets", "mTORC1 signaling", "KRAS signaling up", "Wnt / beta-catenin signaling",
  "Hedgehog signaling", "Notch signaling", "TGF-beta signaling", "Epithelial-mesenchymal transition",
  "Angiogenesis", "Coagulation", "Complement", "Allograft rejection",
  "Estrogen response early", "Androgen response", "Heme metabolism", "Peroxisome",
  "Bile acid metabolism", "Xenobiotic metabolism", "Adipogenesis", "Myogenesis",
  "Spermatogenesis", "Pancreas beta cells", "Apical junction", "UV response up"
)

# A simulated differential-expression experiment.
# Active pathways shift their member genes' test statistics; `frac` = share of members affected.
sim_experiment <- function(n_genes = 10000, seed = 1, scale = 1) {
  set.seed(seed)
  genes <- sprintf("GENE%05d", seq_len(n_genes))
  sizes <- round(exp(runif(length(pathway_names), log(20), log(200))))
  sets <- lapply(sizes, function(s) sample(genes, s))
  names(sets) <- pathway_names
  active <- data.frame(
    set = c("Interferon alpha response", "E2F targets (cell cycle)", "Oxidative phosphorylation", "Hypoxia",
            "Epithelial-mesenchymal transition"),
    shift = c(5, 1.0, -0.8, 0.7, -4.5),          # effect on each affected gene's t statistic
    frac = c(0.5, 0.9, 0.9, 0.8, 0.35)             # share of the set's genes affected
  )
  effect <- setNames(numeric(n_genes), genes)
  # plus a background of strongly changed genes that belong to no particular pathway
  bg <- sample(n_genes, round(0.03 * n_genes))
  effect[bg] <- sample(c(-4.5, 4.5), length(bg), replace = TRUE) * scale
  for (i in seq_len(nrow(active))) {
    m <- sets[[active$set[i]]]
    hit <- sample(m, round(active$frac[i] * length(m)))
    effect[hit] <- effect[hit] + active$shift[i] * scale
  }
  t <- effect + rnorm(n_genes)
  p <- 2 * pnorm(-abs(t))
  data.frame(gene = genes, t = t, logFC = 0.35 * t + rnorm(n_genes, 0, 0.05), p = p,
             padj = p.adjust(p, "BH"), true_effect = effect, stringsAsFactors = FALSE) |>
    structure(sets = sets, active = active)
}

# One-sided hypergeometric (= Fisher's exact) over-representation test.
ora_test <- function(k, K, n, N) {
  phyper(k - 1, K, N - K, n, lower.tail = FALSE)
}

ora_all <- function(de_genes, universe, sets) {
  N <- length(universe); n <- length(intersect(de_genes, universe))
  res <- do.call(rbind, lapply(names(sets), function(s) {
    m <- intersect(sets[[s]], universe); K <- length(m); k <- length(intersect(m, de_genes))
    data.frame(set = s, K = K, k = k, expected = K * n / N, fold = if (K) (k / n) / (K / N) else NA,
               p = if (K) ora_test(k, K, n, N) else NA)
  }))
  res$q <- p.adjust(res$p, "BH")
  res[order(res$p), ]
}

# GSEA running-sum statistic -----------------------------------------------------
# metric: ranking statistic sorted in decreasing order; hits: logical vector, same order.

running_sum <- function(metric, hits, weight = 1) {
  w <- abs(metric)^weight
  NR <- sum(w[hits]); N <- length(metric); NH <- sum(hits)
  step <- ifelse(hits, w / NR, -1 / (N - NH))
  rs <- cumsum(step)
  i <- which.max(abs(rs))
  list(rs = rs, es = rs[i], peak = i)
}

# Fast ES from hit positions only (used for permutations).
es_fast <- function(abs_w, pos, N) {
  pos <- sort(pos); NH <- length(pos)
  ph <- cumsum(abs_w[pos]) / sum(abs_w[pos])
  miss_before <- (pos - seq_len(NH)) / (N - NH)
  after <- ph - miss_before                       # value just after each hit
  before <- c(0, ph[-NH]) - miss_before           # value just before each hit
  mx <- max(after); mn <- min(before)
  if (abs(mx) >= abs(mn)) mx else mn
}

gsea_one <- function(metric_sorted, hits, n_perm = 1000, weight = 1, seed = 1) {
  N <- length(metric_sorted); NH <- sum(hits)
  w <- abs(metric_sorted)^weight
  es <- es_fast(w, which(hits), N)
  set.seed(seed)
  null <- vapply(seq_len(n_perm), function(i) es_fast(w, sample.int(N, NH), N), 0)
  same <- if (es >= 0) null[null >= 0] else null[null < 0]
  nes <- es / mean(abs(same))
  p <- (sum(abs(same) >= abs(es)) + 1) / (length(same) + 1)
  list(es = es, nes = nes, p = p, null = null)
}

gsea_all <- function(stat, sets, n_perm = 1000, weight = 1, seed = 1) {
  o <- order(stat$t, decreasing = TRUE)
  metric <- stat$t[o]; genes <- stat$gene[o]; N <- length(genes)
  w <- abs(metric)^weight
  set.seed(seed)
  # one shared pool of random hit positions per set size keeps this fast
  res <- do.call(rbind, lapply(names(sets), function(s) {
    pos <- match(sets[[s]], genes); pos <- pos[!is.na(pos)]; NH <- length(pos)
    es <- es_fast(w, pos, N)
    null <- vapply(seq_len(n_perm), function(i) es_fast(w, sample.int(N, NH), N), 0)
    same <- if (es >= 0) null[null >= 0] else null[null < 0]
    rs <- running_sum(metric, seq_len(N) %in% pos, weight)
    lead <- if (es >= 0) sum(pos <= rs$peak) else sum(pos >= rs$peak)
    data.frame(set = s, size = NH, es = es, nes = es / mean(abs(same)),
               p = (sum(abs(same) >= abs(es)) + 1) / (length(same) + 1), leading_edge = lead)
  }))
  res$q <- p.adjust(res$p, "BH")
  res[order(res$p, -abs(res$nes)), ]
}
