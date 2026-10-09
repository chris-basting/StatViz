# Data generators and small helpers shared by the regression modules ------------

# Simple linear regression data, optionally breaking one assumption.
sim_linear <- function(n, b0, b1, sd, seed, scenario = "ok") {
  set.seed(seed)
  x <- runif(n, 0, 10)
  e <- switch(
    scenario,
    hetero = rnorm(n, sd = sd * (0.15 + x / 5)),
    skew   = sd * (rexp(n) - 1),
    rnorm(n, sd = sd)
  )
  mu <- b0 + b1 * x
  if (scenario == "curve") mu <- mu + 0.25 * (x - 5)^2 - 2
  y <- mu + e
  if (scenario == "outlier") {
    i <- which.min(abs(x - 9))
    y[i] <- y[i] - 8 * sd - 4
  }
  data.frame(x = x, y = y)
}

# Two predictors with a chosen correlation.
sim_mlr <- function(n, b1, b2, rho, sd, seed) {
  set.seed(seed)
  X <- MASS::mvrnorm(n, c(5, 5), matrix(c(4, 4 * rho, 4 * rho, 4), 2))
  y <- 10 + b1 * X[, 1] + b2 * X[, 2] + rnorm(n, sd = sd)
  data.frame(x1 = X[, 1], x2 = X[, 2], y = y)
}

# Binary outcome from a logistic model.
sim_logistic <- function(n, b0, b1, seed) {
  set.seed(seed)
  x <- runif(n, 0, 10)
  p <- plogis(b0 + b1 * x)
  data.frame(x = x, y = rbinom(n, 1, p))
}

# Grouped data: random intercepts and slopes per group.
sim_mixed <- function(groups, per_group, unbalanced, b0, b1, sd_int, sd_slope, sd_res, seed) {
  set.seed(seed)
  sizes <- if (unbalanced) pmax(2, round(exp(runif(groups, log(2), log(per_group * 3))))) else
    rep(per_group, groups)
  u0 <- rnorm(groups, sd = sd_int)
  u1 <- rnorm(groups, sd = sd_slope)
  g <- rep(seq_len(groups), sizes)
  x <- runif(length(g), 0, 10)
  y <- b0 + u0[g] + (b1 + u1[g]) * x + rnorm(length(g), sd = sd_res)
  d <- data.frame(group = factor(sprintf("School %02d", g)), x = x, y = y)
  attr(d, "truth") <- data.frame(group = levels(d$group), u0 = u0, u1 = u1)
  d
}

fmt <- function(x, d = 2) formatC(x, format = "f", digits = d)

fmt_p <- function(p) ifelse(p < 0.001, "< 0.001", formatC(p, format = "f", digits = 3))

# "p = 0.031" or "p < 0.001"
p_lab <- function(p) ifelse(p < 0.001, "p < 0.001", paste("p =", formatC(p, format = "f", digits = 3)))

# Sequential blue ramp (light = low, dark = high) for continuous magnitudes.
seq_blues <- c("#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b")

# Built-in datasets for the "real data" tabs ------------------------------------

mlr_datasets <- list(
  mtcars = list(name = "Motor Trend cars", data = mtcars, y = "mpg", x = c("wt", "hp"),
                blurb = "32 cars. Fuel economy (mpg) against weight (1000 lb), horsepower and more."),
  swiss = list(name = "Swiss fertility (1888)", data = swiss, y = "Fertility",
               x = c("Education", "Agriculture"),
               blurb = "47 Swiss provinces. Fertility rate against education, farming and religion (% of men/population)."),
  Boston = list(name = "Boston housing", data = MASS::Boston, y = "medv", x = c("rm", "lstat"),
                blurb = "506 neighbourhoods. Median home value ($1000s) against rooms per home (rm), % lower-status population (lstat), crime rate, etc.")
)

logit_datasets <- list(
  pima = list(name = "Pima diabetes", data = transform(MASS::Pima.tr, diabetes = as.integer(type == "Yes"), type = NULL),
              y = "diabetes", x = c("glu", "bmi"),
              blurb = "200 Pima women. Diabetes (1 = yes) against blood glucose (glu), BMI, age, number of pregnancies (npreg), etc."),
  mtcars = list(name = "Motor Trend cars", data = mtcars, y = "am", x = c("wt"),
                blurb = "32 cars. Manual transmission (am = 1) against weight, horsepower, etc."),
  birthwt = list(name = "Low birth weight", data = MASS::birthwt[, c("low", "age", "lwt", "smoke", "ptl", "ht", "ui")],
                 y = "low", x = c("lwt", "smoke"),
                 blurb = "189 births. Low birth weight (1 = under 2.5 kg) against mother's weight in lb (lwt), smoking (smoke), age, etc.")
)
