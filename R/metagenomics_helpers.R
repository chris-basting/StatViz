# Small, readable sequence tools for the metagenomics module ---------------------
# (deliberately simple: forward strand only, toy-sized inputs)

BASES <- c("A", "C", "G", "T")
base_cols <- c(A = "#1baf7a", C = "#2a78d6", G = "#eda100", T = "#e34948")

rand_dna <- function(n, gc = 0.5) {
  p <- c((1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2)
  paste(sample(BASES, n, replace = TRUE, prob = p), collapse = "")
}

chars <- function(s) strsplit(s, "", fixed = TRUE)[[1]]

revcomp <- function(s) paste(rev(chartr("ACGT", "TGCA", chars(s))), collapse = "")

gc_content <- function(s) { x <- chars(s); mean(x %in% c("G", "C")) }

# Introduce random substitutions at rate `rate`.
mutate_dna <- function(s, rate) {
  x <- chars(s)
  hit <- which(runif(length(x)) < rate)
  for (i in hit) x[i] <- sample(setdiff(BASES, x[i]), 1)
  paste(x, collapse = "")
}

# Standard genetic code, codons ordered TTT, TTC, TTA, TTG, TCT, ...
.aa <- strsplit("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "")[[1]]
.cod <- with(expand.grid(third = c("T", "C", "A", "G"), second = c("T", "C", "A", "G"),
                         first = c("T", "C", "A", "G")), paste0(first, second, third))
CODON_TABLE <- setNames(.aa, .cod)

translate <- function(s) {
  n <- nchar(s) %/% 3
  if (n == 0) return("")
  cod <- substring(s, seq(1, by = 3, length.out = n), seq(3, by = 3, length.out = n))
  paste(CODON_TABLE[cod], collapse = "")
}

# Open reading frames: ATG ... stop, on the three forward frames.
find_orfs <- function(s, min_codons = 30) {
  out <- list()
  L <- nchar(s)
  for (f in 0:2) {
    n <- (L - f) %/% 3
    if (n < 2) next
    starts <- seq(1 + f, by = 3, length.out = n)
    cod <- substring(s, starts, starts + 2)
    i <- 1
    while (i <= n) {
      if (cod[i] == "ATG") {
        j <- i
        while (j <= n && !(cod[j] %in% c("TAA", "TAG", "TGA"))) j <- j + 1
        if (j <= n && (j - i) >= min_codons) {
          out[[length(out) + 1]] <- data.frame(frame = f + 1, start = starts[i], end = starts[j] + 2,
                                               codons = j - i)
        }
        i <- j + 1
      } else i <- i + 1
    }
  }
  if (!length(out)) return(data.frame(frame = integer(0), start = integer(0), end = integer(0), codons = integer(0)))
  do.call(rbind, out)
}

# Random protein-coding gene of `codons` codons (ATG ... no internal stops ... stop).
rand_gene <- function(codons, gc = 0.5) {
  sense <- names(CODON_TABLE)[CODON_TABLE != "*" & names(CODON_TABLE) != "ATG"]
  w <- vapply(sense, function(c) prod(ifelse(chars(c) %in% c("G", "C"), gc, 1 - gc)), 0)
  paste0("ATG", paste(sample(sense, codons - 1, replace = TRUE, prob = w), collapse = ""), "TAA")
}

# Pairwise alignment by dynamic programming ------------------------------------
# type = "global" (Needleman-Wunsch) or "local" (Smith-Waterman).
dp_align <- function(a, b, match = 2, mismatch = -1, gap = -2, type = "global") {
  x <- chars(a); y <- chars(b); n <- length(x); m <- length(y)
  S <- matrix(0, n + 1, m + 1); P <- matrix("", n + 1, m + 1)   # P = pointer: d(iag), u(p), l(eft)
  if (type == "global") {
    S[, 1] <- gap * (0:n); S[1, ] <- gap * (0:m)
    P[-1, 1] <- "u"; P[1, -1] <- "l"
  }
  for (i in seq_len(n)) for (j in seq_len(m)) {
    s <- if (x[i] == y[j]) match else mismatch
    cand <- c(d = S[i, j] + s, u = S[i, j + 1] + gap, l = S[i + 1, j] + gap)
    if (type == "local") cand <- c(cand, z = 0)
    k <- which.max(cand)
    S[i + 1, j + 1] <- cand[k]; P[i + 1, j + 1] <- names(cand)[k]
  }
  # traceback
  if (type == "global") { i <- n + 1; j <- m + 1 } else {
    w <- which(S == max(S), arr.ind = TRUE)[1, ]; i <- w[1]; j <- w[2]
  }
  path <- list(); ax <- character(0); ay <- character(0)
  repeat {
    path[[length(path) + 1]] <- c(i, j)
    if (type == "local" && (S[i, j] == 0 || P[i, j] %in% c("", "z"))) break
    if (type == "global" && i == 1 && j == 1) break
    p <- P[i, j]
    if (p == "d") { ax <- c(x[i - 1], ax); ay <- c(y[j - 1], ay); i <- i - 1; j <- j - 1 }
    else if (p == "u") { ax <- c(x[i - 1], ax); ay <- c("-", ay); i <- i - 1 }
    else if (p == "l") { ax <- c("-", ax); ay <- c(y[j - 1], ay); j <- j - 1 }
    else break
  }
  path <- do.call(rbind, rev(path))
  mid <- ifelse(ax == ay, "|", ifelse(ax == "-" | ay == "-", " ", "."))
  list(S = S, path = path, a = paste(ax, collapse = ""), b = paste(ay, collapse = ""),
       mid = paste(mid, collapse = ""), score = if (type == "global") S[n + 1, m + 1] else max(S),
       identity = if (length(ax)) mean(ax == ay) else 0, x = x, y = y)
}

# Fast local alignment score for proteins (identity scoring), used for the toy homology search.
local_score <- function(a, b, match = 3, mismatch = -1, gap = -3) {
  x <- chars(a); y <- chars(b); m <- length(y)
  prev <- numeric(m + 1); best <- 0
  for (i in seq_along(x)) {
    cur <- numeric(m + 1)
    s <- ifelse(y == x[i], match, mismatch)
    for (j in seq_len(m)) {
      v <- max(0, prev[j] + s[j], prev[j + 1] + gap, cur[j] + gap)
      cur[j + 1] <- v
    }
    best <- max(best, cur); prev <- cur
  }
  best
}

# de Bruijn assembly ------------------------------------------------------------

# Reads sampled uniformly from the forward strand, with substitution errors.
sample_reads <- function(genome, read_len, coverage, error = 0) {
  L <- nchar(genome)
  n <- max(1, round(coverage * L / read_len))
  starts <- sample.int(L - read_len + 1, n, replace = TRUE)
  reads <- substring(genome, starts, starts + read_len - 1)
  if (error > 0) reads <- vapply(reads, mutate_dna, "", rate = error, USE.NAMES = FALSE)
  data.frame(start = starts, read = reads)
}

count_kmers <- function(reads, k) {
  km <- unlist(lapply(reads, function(r) {
    L <- nchar(r); if (L < k) return(character(0))
    substring(r, 1:(L - k + 1), k:L)
  }))
  sort(table(km), decreasing = TRUE)
}

# Build the de Bruijn graph from k-mers and compact it into unitigs (contigs).
debruijn_assemble <- function(kmers) {
  k <- nchar(kmers[1])
  pre <- substr(kmers, 1, k - 1); suf <- substr(kmers, 2, k)
  out_e <- split(seq_along(kmers), pre); in_e <- split(seq_along(kmers), suf)
  nodes <- unique(c(pre, suf))
  outdeg <- setNames(lengths(out_e[nodes]), nodes); outdeg[is.na(outdeg)] <- 0
  indeg <- setNames(lengths(in_e[nodes]), nodes); indeg[is.na(indeg)] <- 0
  simple <- function(v) indeg[[v]] == 1 && outdeg[[v]] == 1
  used <- logical(length(kmers)); unitigs <- character(0); u_start <- character(0); u_end <- character(0)

  walk <- function(e) {
    seq <- kmers[e]; used[e] <<- TRUE; v <- suf[e]
    while (simple(v)) {
      e2 <- out_e[[v]]
      if (used[e2]) break
      used[e2] <<- TRUE; seq <- paste0(seq, substr(kmers[e2], k, k)); v <- suf[e2]
    }
    list(seq = seq, end = v)
  }
  for (v in nodes[!vapply(nodes, simple, TRUE)]) for (e in out_e[[v]]) if (!used[e]) {
    w <- walk(e); unitigs <- c(unitigs, w$seq); u_start <- c(u_start, v); u_end <- c(u_end, w$end)
  }
  while (any(!used)) {                       # isolated cycles
    e <- which(!used)[1]; w <- walk(e)
    unitigs <- c(unitigs, w$seq); u_start <- c(u_start, pre[e]); u_end <- c(u_end, w$end)
  }
  # unitig graph: A -> B when A ends where B starts
  edges <- do.call(rbind, lapply(seq_along(unitigs), function(a) {
    b <- which(u_start == u_end[a]); if (length(b)) data.frame(from = a, to = b) else NULL
  }))
  list(contigs = unitigs, k = k, n_nodes = length(nodes), n_kmers = length(kmers),
       edges = if (is.null(edges)) data.frame(from = integer(0), to = integer(0)) else edges,
       branching = sum(indeg > 1 | outdeg > 1))
}

n50 <- function(lengths) {
  l <- sort(lengths, decreasing = TRUE)
  l[which(cumsum(l) >= sum(l) / 2)[1]]
}

# Simulated contigs for binning: GC content and coverage per genome.
sim_contigs <- function(n_genomes, similarity, mean_len, seed) {
  set.seed(seed)
  gc <- c(0.35, 0.62, 0.48, 0.55)[seq_len(n_genomes)]
  cov <- c(40, 12, 90, 25)[seq_len(n_genomes)]
  # similarity pulls genome GC and coverage toward a common value
  gc <- gc + similarity * (0.5 - gc)
  cov <- exp(log(cov) + similarity * (log(30) - log(cov)))
  do.call(rbind, lapply(seq_len(n_genomes), function(g) {
    n <- 40
    len <- round(rlnorm(n, log(mean_len), 0.6))
    data.frame(genome = paste("Genome", LETTERS[g]), length = len,
               gc = rbinom(n, len, gc[g]) / len,
               coverage = rgamma(n, shape = 25, rate = 25 / cov[g]) * rpois(n, 1000) / 1000)
  }))
}

taxa <- data.frame(
  taxon = c("Bacteroides fragilis", "Escherichia coli", "Lactobacillus reuteri", "Prevotella copri",
            "Unknown species"),
  genome_mb = c(5.2, 5.0, 2.0, 3.5, 4.0),
  in_db = c(TRUE, TRUE, TRUE, TRUE, FALSE)
)
