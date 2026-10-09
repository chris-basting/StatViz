# Shared pieces for the hypothesis-test modules (t-test, ANOVA, Wilcoxon) -------

two_cols <- c(Control = "#2a78d6", Treatment = "#eb6834")

# Random noise with a chosen shape, scaled to roughly the requested SD.
draw_noise <- function(n, sd, shape = "normal") {
  switch(
    shape,
    skewed = sd * (rexp(n) - 1),                  # long right tail
    heavy  = sd * rt(n, df = 2) / 1.5,            # occasional extreme values
    rnorm(n, sd = sd)
  )
}

shape_choices <- c("Bell-shaped (normal)" = "normal", "Skewed (long right tail)" = "skewed",
                   "Heavy tails (occasional extreme values)" = "heavy")

# Two independent groups.
sim_two_groups <- function(n, diff, sd, seed, shape = "normal", base = 10) {
  set.seed(seed)
  data.frame(
    group = factor(rep(c("Control", "Treatment"), each = n), levels = c("Control", "Treatment")),
    y = c(base + draw_noise(n, sd, shape), base + diff + draw_noise(n, sd, shape))
  )
}

# k groups whose true means are spread out by `spread`.
sim_k_groups <- function(k, n, spread, sd, seed, shape = "normal", base = 20) {
  set.seed(seed)
  mu <- base + spread * seq(-1, 1, length.out = k) * sqrt(3)   # SD of the means ~ spread
  if (k == 1) mu <- base
  g <- factor(rep(LETTERS[seq_len(k)], each = n))
  data.frame(group = g, y = mu[as.integer(g)] + draw_noise(n * k, sd, shape))
}

# Paired before/after measurements on the same people.
sim_paired <- function(n, change, sd_person, sd_noise, seed, shape = "normal") {
  set.seed(seed)
  person <- rnorm(n, 50, sd_person)
  data.frame(id = seq_len(n),
             before = person + draw_noise(n, sd_noise, shape),
             after = person + change + draw_noise(n, sd_noise, shape))
}

# Density curve of a null distribution with the "as or more extreme" region shaded.
# side = "two" shades both tails beyond |obs|; "upper" shades above obs.
null_curve_plot <- function(x, dens, obs, side = "two", xlab = NULL, p = NULL, ymax = NULL) {
  ymax <- ymax %||% max(dens)
  df <- data.frame(x = x, d = dens)
  tail <- if (side == "two") abs(df$x) >= abs(obs) else df$x >= obs
  df$left <- ifelse(tail & df$x < 0, df$d, 0)     # two areas so the shading doesn't bridge the gap
  df$right <- ifelse(tail & df$x >= 0, df$d, 0)
  g <- ggplot(df, aes(x, d)) +
    geom_area(aes(y = left), fill = sv_col$accent, alpha = 0.45) +
    geom_area(aes(y = right), fill = sv_col$accent, alpha = 0.45) +
    geom_line(colour = sv_col$primary, linewidth = 1.1) +
    geom_vline(xintercept = obs, colour = sv_col$ink, linewidth = 0.9)
  if (side == "two") g <- g + geom_vline(xintercept = -obs, colour = sv_col$ink, linewidth = 0.5, linetype = "22")
  if (!is.null(p)) {
    g <- g + annotate("label", x = obs, y = ymax * 0.9, label = paste0("observed\n", p_lab(p)),
                      hjust = if (obs > mean(range(x))) 1.05 else -0.05, colour = sv_col$ink,
                      fill = sv_col$surface, label.size = 0, size = 3.8)
  }
  g + scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    coord_cartesian(ylim = c(0, ymax * 1.02)) +
    labs(x = xlab, y = "Probability density") + theme_statviz()
}

# Points per group with mean and 95% CI, groups on the x axis (position = identity).
group_points_plot <- function(d, y = "y", group = "group", ylab = "Value", colours = NULL, show_ci = TRUE,
                              jitter = 0.12) {
  s <- do.call(rbind, lapply(split(d[[y]], d[[group]]), function(v) {
    m <- mean(v); se <- sd(v) / sqrt(length(v)); h <- qt(0.975, max(1, length(v) - 1)) * se
    data.frame(m = m, lo = m - h, hi = m + h)
  }))
  s$g <- factor(rownames(s), levels = levels(factor(d[[group]])))
  set.seed(3)
  g <- ggplot(d, aes(.data[[group]], .data[[y]]))
  g <- if (is.null(colours)) {
    g + geom_jitter(width = jitter, height = 0, colour = sv_col$primary, alpha = 0.6, size = 2.2)
  } else {
    g + geom_jitter(aes(colour = .data[[group]]), width = jitter, height = 0, alpha = 0.7, size = 2.2) +
      scale_colour_manual(values = colours, guide = "none")
  }
  if (show_ci) {
    g <- g + geom_errorbar(data = s, aes(x = g, ymin = lo, ymax = hi), inherit.aes = FALSE,
                           width = 0.08, colour = sv_col$ink, linewidth = 0.8,
                           position = position_nudge(x = 0.28))
  }
  g + geom_point(data = s, aes(x = g, y = m), inherit.aes = FALSE, colour = sv_col$ink, size = 3.5,
                 shape = 18, position = position_nudge(x = 0.28)) +
    labs(x = NULL, y = ylab) + theme_statviz()
}

# Plain-language "Results" sentence fragments.
p_sentence <- function(p, alpha = 0.05) {
  if (p < alpha) sprintf("%s, below %s: the difference is unlikely to be just chance.", p_lab(p), alpha)
  else sprintf("%s, not below %s: a difference this size could easily be chance.", p_lab(p), alpha)
}

# Built-in datasets ------------------------------------------------------------

tt_datasets <- list(
  sleep = list(name = "Sleep drugs (paired)", data = transform(sleep, group = factor(group, labels = c("Drug 1", "Drug 2"))),
               y = "extra", group = "group", id = "ID",
               blurb = "10 patients each tried two sleep drugs. extra = extra hours of sleep. The same people got both drugs, so the data are paired."),
  ToothGrowth = list(name = "Guinea-pig tooth growth", data = ToothGrowth, y = "len", group = "supp", id = NULL,
                     blurb = "60 guinea pigs given vitamin C as orange juice (OJ) or ascorbic acid (VC). len = tooth length."),
  mtcars = list(name = "Car fuel economy", data = transform(mtcars, am = factor(am, labels = c("Automatic", "Manual"))),
                y = "mpg", group = "am", id = NULL,
                blurb = "32 cars: miles per gallon for automatic vs manual transmissions."),
  PlantGrowth = list(name = "Plant growth", data = PlantGrowth, y = "weight", group = "group", id = NULL,
                     blurb = "30 plants: a control and two treatments. Pick any two groups to compare."),
  chickwts = list(name = "Chick weights by feed", data = chickwts, y = "weight", group = "feed", id = NULL,
                  blurb = "71 chicks raised on six different feeds. Pick any two feeds to compare.")
)

anova_datasets <- list(
  PlantGrowth = list(name = "Plant growth", data = PlantGrowth, y = "weight", a = "group", b = NULL,
                     blurb = "Dried weight of plants under a control and two treatments."),
  chickwts = list(name = "Chick weights by feed", data = chickwts, y = "weight", a = "feed", b = NULL,
                  blurb = "Weight of chicks after six weeks on one of six feeds."),
  InsectSprays = list(name = "Insect sprays", data = InsectSprays, y = "count", a = "spray", b = NULL,
                      blurb = "Number of insects found on plots treated with six sprays."),
  ToothGrowth = list(name = "Tooth growth (two factors)", data = transform(ToothGrowth, dose = factor(dose)),
                     y = "len", a = "dose", b = "supp",
                     blurb = "Tooth length by vitamin C dose (0.5, 1, 2 mg) and delivery method (OJ or VC)."),
  warpbreaks = list(name = "Yarn breaks (two factors)", data = warpbreaks, y = "breaks", a = "tension", b = "wool",
                    blurb = "Number of warp breaks per loom by wool type (A, B) and tension (L, M, H).")
)

wx_datasets <- list(
  mtcars = list(name = "Car fuel economy (2 groups)", type = "two",
                data = transform(mtcars, am = factor(am, labels = c("Automatic", "Manual"))), y = "mpg", group = "am",
                blurb = "Miles per gallon for automatic vs manual cars."),
  ToothGrowth = list(name = "Tooth growth (2 groups)", type = "two", data = ToothGrowth, y = "len", group = "supp",
                     blurb = "Tooth length with orange juice (OJ) vs ascorbic acid (VC)."),
  sleep = list(name = "Sleep drugs (paired)", type = "paired",
               data = transform(sleep, group = factor(group, labels = c("Drug 1", "Drug 2"))),
               y = "extra", group = "group", id = "ID",
               blurb = "10 patients, each tried both drugs. The signed-rank test compares the drugs within each person."),
  InsectSprays = list(name = "Insect sprays (6 groups)", type = "k", data = InsectSprays, y = "count", group = "spray",
                      blurb = "Insect counts under six sprays: counts are skewed, a natural case for Kruskal-Wallis."),
  airquality = list(name = "Ozone by month (5 groups)", type = "k",
                    data = transform(na.omit(airquality[, c("Ozone", "Month")]), Month = factor(month.abb[Month], levels = month.abb[5:9])),
                    y = "Ozone", group = "Month",
                    blurb = "New York ozone readings, May-September 1973. Strongly skewed, with some very high days.")
)
