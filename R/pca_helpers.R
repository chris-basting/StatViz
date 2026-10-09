# PCA computations, written out by hand so the app can show every step --------

# Fit PCA via eigen-decomposition of the covariance (or correlation) matrix.
pca_fit <- function(X, scale = TRUE) {
  X <- as.matrix(X)
  means <- colMeans(X)
  sds <- apply(X, 2, sd)
  Z <- sweep(X, 2, means)
  if (scale) Z <- sweep(Z, 2, sds, "/")
  S <- crossprod(Z) / (nrow(Z) - 1)
  e <- eigen(S, symmetric = TRUE)
  V <- e$vectors
  # Eigenvectors are only defined up to sign; make the largest loading positive
  # so the picture doesn't flip around between datasets.
  for (j in seq_len(ncol(V))) {
    if (V[which.max(abs(V[, j])), j] < 0) V[, j] <- -V[, j]
  }
  pcs <- paste0("PC", seq_len(ncol(V)))
  dimnames(V) <- list(colnames(X), pcs)
  scores <- Z %*% V
  colnames(scores) <- pcs
  values <- pmax(e$values, 0)
  list(X = X, Z = Z, means = means, sds = sds, scale = scale, S = S,
       values = setNames(values, pcs), V = V, scores = scores,
       prop = setNames(values / sum(values), pcs))
}

# Rebuild the original data from the first k components.
pca_reconstruct <- function(fit, k) {
  Zhat <- fit$scores[, seq_len(k), drop = FALSE] %*% t(fit$V[, seq_len(k), drop = FALSE])
  if (fit$scale) Zhat <- sweep(Zhat, 2, fit$sds, "*")
  sweep(Zhat, 2, fit$means, "+")
}

# Correlated 2-D cloud for the intuition tab.
gen_2d <- function(n, rho, sd_x, sd_y, seed) {
  set.seed(seed)
  Sigma <- matrix(c(sd_x^2, rho * sd_x * sd_y, rho * sd_x * sd_y, sd_y^2), 2)
  X <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = Sigma)
  colnames(X) <- c("x1", "x2")
  X
}

# A 3-D cloud that lies near a (randomly tilted) plane.
gen_3d <- function(n, thickness, seed) {
  set.seed(seed)
  latent <- cbind(rnorm(n, sd = 3), rnorm(n, sd = 1.5), rnorm(n, sd = thickness))
  R <- qr.Q(qr(matrix(rnorm(9), 3)))
  X <- latent %*% t(R)
  colnames(X) <- c("x1", "x2", "x3")
  X
}

# Variance of the data projected onto unit vector at `deg` degrees.
proj_variance <- function(S, deg) {
  th <- deg * pi / 180
  u <- c(cos(th), sin(th))
  drop(t(u) %*% S %*% u)
}

# Built-in datasets ------------------------------------------------------------

pca_datasets <- list(
  iris = list(
    name = "Iris flowers",
    data = iris, group = "Species",
    blurb = "150 flowers, 4 measurements (cm). Do the three species separate along the first PC?"
  ),
  USArrests = list(
    name = "US arrests by state",
    data = USArrests, group = NULL,
    blurb = "Crime rates per 100,000 and % urban population for 50 states. A classic demo of why scaling matters: try unchecking 'Scale'."
  ),
  mtcars = list(
    name = "Motor Trend cars",
    data = transform(mtcars, cyl = factor(cyl), vs = NULL, am = NULL),
    group = "cyl",
    blurb = "32 cars, 9 numeric specs. Colored by number of cylinders."
  ),
  swiss = list(
    name = "Swiss fertility (1888)",
    data = swiss, group = NULL,
    blurb = "47 French-speaking Swiss provinces: fertility and socio-economic indicators."
  )
)

pca_dataset_choices <- c(
  setNames(names(pca_datasets), vapply(pca_datasets, `[[`, "", "name")),
  "Upload a CSV..." = "upload"
)

# Long-format matrix for heatmaps.
mat_long <- function(M) {
  data.frame(
    row = factor(rep(rownames(M), times = ncol(M)), levels = rev(rownames(M))),
    col = factor(rep(colnames(M), each = nrow(M)), levels = colnames(M)),
    value = as.vector(M)
  )
}

mat_heatmap <- function(M, fill_name = NULL, digits = 2) {
  df <- mat_long(M)
  lim <- max(abs(df$value), 1e-9)
  ggplot(df, aes(col, row, fill = value)) +
    geom_tile(colour = sv_col$surface, linewidth = 1.5) +
    geom_text(aes(label = formatC(value, format = "f", digits = digits)),
              colour = sv_col$ink, size = 3.6) +
    scale_fill_diverging(lim, fill_name) +
    scale_x_discrete(position = "top") +
    labs(x = NULL, y = NULL) +
    theme_statviz() +
    theme(panel.grid = element_blank(), legend.position = "right",
          axis.text.x = element_text(angle = 30, hjust = 0))
}
