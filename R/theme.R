# Shared colors and ggplot theme for every StatViz module ---------------------

sv_col <- list(
  ink       = "#0b0b0b",
  ink2      = "#52514e",
  muted     = "#8a8984",
  grid      = "#e6e5e1",
  neutral   = "#f0efec",
  surface   = "#fcfcfb",
  primary   = "#2a78d6",  # blue  - main series / "the thing we're learning"
  accent    = "#eb6834",  # orange - the learner's guess / projections
  good      = "#1baf7a",  # aqua  - the optimal answer
  neg       = "#2a78d6",  # diverging: negative pole
  pos       = "#e34948"   # diverging: positive pole
)

# Fixed categorical order (validated). Scatterplots are capped at 3 groups for
# colorblind safety; beyond that we fall back to the full 8 with a warning.
sv_groups <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100",
               "#e87ba4", "#008300", "#4a3aa7", "#e34948")

theme_statviz <- function(base_size = 13) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.background  = element_rect(fill = sv_col$surface, colour = NA),
      panel.background = element_rect(fill = sv_col$surface, colour = NA),
      panel.grid.major = element_line(colour = sv_col$grid, linewidth = 0.4),
      panel.grid.minor = element_blank(),
      axis.title       = element_text(colour = sv_col$ink2),
      axis.text        = element_text(colour = sv_col$ink2),
      plot.title       = element_text(colour = sv_col$ink, face = "bold"),
      plot.subtitle    = element_text(colour = sv_col$ink2),
      legend.position  = "top",
      legend.title     = element_text(colour = sv_col$ink2),
      legend.text      = element_text(colour = sv_col$ink)
    )
}

scale_fill_diverging <- function(limit, name = NULL) {
  scale_fill_gradient2(low = sv_col$neg, mid = sv_col$neutral, high = sv_col$pos,
                       midpoint = 0, limits = c(-limit, limit), name = name)
}

group_palette <- function(levels) {
  setNames(sv_groups[seq_along(levels)], levels)
}
