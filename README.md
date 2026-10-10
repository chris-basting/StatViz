# StatViz

Live app: https://chris-basting.shinyapps.io/statviz/

An interactive R Shiny app for learning statistical and data-science methods by
moving sliders and watching the plots react.

## Run

```r
# install.packages(c("shiny", "bslib", "ggplot2", "plotly", "ggrepel", "MASS", "lme4", "igraph"))
shiny::runApp()
```

## Modules

Every page explains each idea in plain language first, then shows the formula with an
"In words" translation; heavier derivations sit in collapsible sections.

### Principal Component Analysis
The big idea (rotate a line), step by step (center → covariance → eigenvectors → scores),
how many PCs, interpreting loadings, reconstruction, 3-D, quiz.

### Regression (menu)
- **Simple linear regression**: fit a line by hand with an RSS landscape, slope/intercept
  formulas and R², sampling variability and confidence/prediction bands, assumption checks
  (LINE), click-to-add outliers with leverage and Cook's distance.
- **Multiple linear regression**: the fitted plane in 3-D, "holding others constant"
  (confounding, added-variable plot), multicollinearity and VIF, real data (mtcars, swiss,
  Boston) with plain-language coefficient sentences.
- **Logistic regression**: why not a straight line, odds and log-odds, maximum likelihood
  by hand with a likelihood landscape, thresholds/confusion matrix/ROC, real data (Pima
  diabetes, mtcars, birth weight) with odds ratios.
- **Mixed-effects regression**: complete vs no vs partial pooling, shrinkage, ICC and
  variance components, a coverage simulation showing why ignoring groups is overconfident,
  and the lme4 sleep study.

### Comparing groups (menu)
- **t-tests**: signal vs noise, where p-values come from (null distribution + simulation),
  paired vs unpaired, Type I/II errors and power, real data (sleep, ToothGrowth, mtcars,
  PlantGrowth, chickwts).
- **ANOVA**: why not many t-tests (false alarms pile up), between vs within sums of squares
  and the F ratio, the F distribution, Tukey post-hoc comparisons, two-way ANOVA with an
  interaction plot, real data.
- **Wilcoxon & rank tests**: ranks instead of values (and the effect of an outlier), exact
  and shuffled null distributions, t-test vs Wilcoxon power under different noise shapes,
  signed-rank for paired data, Kruskal-Wallis, real data with both tests side by side.

### Metagenomics
The pipeline and shotgun coverage (Lander-Waterman), read quality and trimming (Phred scores),
pairwise alignment by dynamic programming (global and local, with the full score table and
traceback), de Bruijn graph assembly from simulated reads (k-mers, repeats, sequencing errors,
N50), binning contigs into MAGs by GC content and coverage, annotation (ORF finding plus a toy
homology search), and taxonomic profiling (k-mer classification, genome-size bias, diversity).
The assembler, aligner and gene finder are small working implementations in base R.

## Layout

```
app.R              # navbar shell; one nav_panel per method
R/theme.R              # shared palette + ggplot theme
R/ui_helpers.R         # explain(), math_words(), deeper(), axes_note(), quiz helpers
R/pca_helpers.R        # PCA math (by hand) + data generators
R/regression_helpers.R # simulators + built-in datasets for the regression pages
R/test_helpers.R       # simulators, datasets and null-distribution plots for the test pages
R/metagenomics_helpers.R # toy sequence tools: translation, ORFs, alignment, de Bruijn assembly
R/mod_pca.R            # PCA module
R/mod_linreg.R         # simple linear regression module
R/mod_mlr.R            # multiple linear regression module
R/mod_logreg.R         # logistic regression module
R/mod_mixed.R          # mixed-effects regression module
R/mod_ttest.R          # t-test module
R/mod_anova.R          # ANOVA module
R/mod_wilcox.R         # Wilcoxon / rank tests module
R/mod_metagenomics.R   # metagenomics module
www/styles.css
www/rotate.js      # client-side Rotate/Pause animation
```

To add a method, create `R/mod_<method>.R` with `<method>_ui()` / `<method>_server()`
and add it to `pages` and the navbar in `app.R`.
