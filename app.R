# StatViz: interactive explanations of statistics & data-science methods
#
# Run with:  shiny::runApp()
# Files in R/ are sourced automatically by Shiny.

library(shiny)
library(bslib)
library(ggplot2)
library(plotly)

# Every method page: nav value, menu label, icon, blurb for the home page.
pages <- list(
  list(id = "pca", label = "Principal Component Analysis", icon = "compress",
       blurb = "Find the directions that capture the most variation, and use them to compress and visualize many variables at once."),
  list(id = "linreg", label = "Simple linear regression", icon = "chart-line",
       blurb = "Fit a line by hand, see what least squares minimizes, and learn how sure you can be about the slope."),
  list(id = "mlr", label = "Multiple linear regression", icon = "cube",
       blurb = "Several predictors at once: what 'holding the others constant' really means, confounding and multicollinearity."),
  list(id = "logreg", label = "Logistic regression", icon = "toggle-on",
       blurb = "Predict yes/no outcomes: probabilities, odds and log-odds, maximum likelihood, and turning probabilities into decisions."),
  list(id = "mixed", label = "Mixed-effects regression", icon = "layer-group",
       blurb = "Grouped data like students in schools: random effects, shrinkage, and why ignoring groups makes you overconfident."),
  list(id = "ttest", label = "t-tests", icon = "scale-balanced",
       blurb = "Is a difference between two groups real? Signal vs noise, where p-values come from, paired designs and power."),
  list(id = "anova", label = "ANOVA", icon = "chart-column",
       blurb = "Compare several groups at once: between vs within variation, the F ratio, post-hoc tests and interactions."),
  list(id = "wilcox", label = "Wilcoxon & rank tests", icon = "arrow-up-1-9",
       blurb = "Tests built on ranks: rank-sum, signed-rank and Kruskal-Wallis, and when they beat the t-test.")
)

home_card <- function(pg) {
  card(
    card_header(icon(pg$icon), pg$label),
    p(pg$blurb),
    actionButton(paste0("go_", pg$id), "Start", class = "btn-primary", width = "120px")
  )
}

ui <- page_navbar(
  title = tagList(icon("chart-line"), "StatViz"),
  id = "main_nav",
  theme = bs_theme(
    version = 5,
    primary = sv_col$primary,
    bg = "#ffffff", fg = sv_col$ink,
    base_font = font_collection("Inter", "system-ui", "-apple-system", "Segoe UI", "sans-serif")
  ),
  header = tags$head(
    tags$link(rel = "stylesheet", href = "styles.css"),
    # Close the Regression menu after jumping to one of its pages from Home.
    tags$script(HTML("
      $(document).on('shown.bs.tab', '.navbar', function () {
        setTimeout(function () {
          $('.navbar .dropdown-menu.show').removeClass('show');
          $('.navbar .dropdown-toggle.show').removeClass('show').attr('aria-expanded', 'false');
        }, 50);
      });"))
  ),

  nav_panel(
    "Home", value = "home", icon = icon("house"),
    div(
      class = "container py-4", style = "max-width: 1000px;",
      h1("StatViz"),
      p(class = "lead",
        "Hands-on explanations of statistical and data-science methods.",
        "Every idea comes with sliders to move and plots that react, so you can build intuition",
        "before (or alongside) the math."),
      h4(class = "mt-4", "Dimension reduction"),
      layout_column_wrap(width = 1 / 2, home_card(pages[[1]])),
      h4(class = "mt-4", "Regression"),
      p(class = "text-muted", "Suggested order: simple → multiple → logistic → mixed effects."),
      do.call(layout_column_wrap, c(list(width = 1 / 2), lapply(pages[2:5], home_card))),
      h4(class = "mt-4", "Comparing groups"),
      p(class = "text-muted", "Suggested order: t-tests → ANOVA → Wilcoxon & rank tests."),
      do.call(layout_column_wrap, c(list(width = 1 / 2), lapply(pages[6:8], home_card)))
    )
  ),

  nav_panel("PCA", value = "pca", icon = icon("compress"), pca_ui("pca")),

  nav_menu(
    "Regression", icon = icon("chart-line"),
    nav_panel("Simple linear regression", value = "linreg", linreg_ui("linreg")),
    nav_panel("Multiple linear regression", value = "mlr", mlr_ui("mlr")),
    nav_panel("Logistic regression", value = "logreg", logreg_ui("logreg")),
    nav_panel("Mixed-effects regression", value = "mixed", mixed_ui("mixed"))
  ),

  nav_menu(
    "Comparing groups", icon = icon("scale-balanced"),
    nav_panel("t-tests", value = "ttest", ttest_ui("ttest")),
    nav_panel("ANOVA", value = "anova", anova_ui("anova")),
    nav_panel("Wilcoxon & rank tests", value = "wilcox", wilcox_ui("wilcox"))
  )
)

server <- function(input, output, session) {
  lapply(pages, function(pg) observeEvent(input[[paste0("go_", pg$id)]], nav_select("main_nav", pg$id)))
  pca_server("pca")
  linreg_server("linreg")
  mlr_server("mlr")
  logreg_server("logreg")
  mixed_server("mixed")
  ttest_server("ttest")
  anova_server("anova")
  wilcox_server("wilcox")
}

shinyApp(ui, server)
