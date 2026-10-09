# Simple linear regression module ----------------------------------------------
#
# Tabs:
#   1. The best line      - fit a line by hand; residuals, squares, RSS landscape
#   2. Formulas & R²      - slope/intercept formulas; SST = SSR + SSE
#   3. How sure are we?   - sampling variability, SE, CI, prediction bands
#   4. Checking assumptions - residual plots under broken assumptions
#   5. Outliers & influence - click to add points
#   6. Quiz

linreg_ui <- function(id) {
  ns <- NS(id)

  layout_sidebar(
    sidebar = sidebar(
      width = 300,
      show_on(ns, c("fit", "formula", "uncert"),
        h6("Simulated data"),
        p(class = "small text-muted",
          "We make up data from a known 'true' line plus random noise, so we can check how well",
          "regression recovers the truth."),
        sliderInput(ns("n"), "Number of points", 5, 200, 30, step = 5),
        sliderInput(ns("b1_true"), "True slope", -2, 2, 0.8, step = 0.1),
        sliderInput(ns("sd"), "Noise (SD of the scatter around the line)", 0.2, 5, 1.5, step = 0.1),
        actionButton(ns("resample"), "Draw a new sample", icon = icon("shuffle"),
                     class = "btn-outline-secondary btn-sm w-100")
      ),
      show_on(ns, "fit",
        hr(),
        h6("Your line"),
        sliderInput(ns("my_b0"), "Intercept (where the line crosses x = 0)", -5, 12, 5, step = 0.1),
        sliderInput(ns("my_b1"), "Slope (change in y per 1 unit of x)", -2, 3, 0, step = 0.05),
        checkboxInput(ns("show_resid"), "Show residuals", TRUE),
        checkboxInput(ns("show_sq"), "Show squared residuals", FALSE),
        actionButton(ns("reveal"), "Show the least-squares line", icon = icon("bullseye"),
                     class = "btn-primary w-100")
      ),
      show_on(ns, "uncert",
        hr(),
        checkboxInput(ns("show_ci"), "Confidence band (where the true line probably is)", TRUE),
        checkboxInput(ns("show_pi"), "Prediction band (where new points probably land)", FALSE)
      ),
      show_on(ns, "assume",
        selectInput(ns("scenario"), "Scenario", c(
          "Everything fine" = "ok",
          "Curved relationship" = "curve",
          "Spread grows with x" = "hetero",
          "Skewed errors" = "skew",
          "One outlier" = "outlier"
        )),
        sliderInput(ns("n_a"), "Number of points", 20, 300, 100, step = 10),
        actionButton(ns("resample_a"), "Draw a new sample", icon = icon("shuffle"),
                     class = "btn-outline-secondary btn-sm w-100"),
        uiOutput(ns("assume_advice"))
      ),
      show_on(ns, "influence",
        p("Click anywhere on the plot to add a point. Watch how much (or how little) the line moves."),
        actionButton(ns("undo"), "Remove last point", icon = icon("rotate-left"),
                     class = "btn-outline-secondary btn-sm w-100 mb-2"),
        actionButton(ns("reset"), "Remove all added points", icon = icon("trash"),
                     class = "btn-outline-secondary btn-sm w-100")
      ),
      show_on(ns, "quiz",
        p("Test yourself on the ideas from the other tabs.")
      )
    ),

    navset_underline(
      id = ns("tab"),

      # ==== 1. The best line ===================================================
      nav_panel(
        "1 · The best line", value = "fit",
        explain(withMathJax(
          p(strong("The goal:"), "use one variable (x) to predict another (y) with a straight line.",
            "A line is described by two numbers: its ", strong("intercept"),
            " (the predicted y when x is 0) and its ", strong("slope"),
            " (how much the prediction goes up for each 1-unit step in x)."),
          math_words(
            "$$\\hat{y} = b_0 + b_1\\,x$$",
            "the predicted y equals the intercept plus the slope times x.",
            symbols = list("\\(\\hat{y}\\)" = "predicted y (the 'hat' means 'estimate')",
                           "\\(b_0\\)" = "intercept", "\\(b_1\\)" = "slope")
          ),
          p(strong("Try it:"), "move the intercept and slope sliders to fit the line by eye.",
            "For each point, the vertical gap between the point and the line is its ", strong("residual"),
            ": actual y minus predicted y. (Unlike PCA, the gaps are vertical, because we're",
            "predicting y and only care about mistakes in y.)"),
          math_words(
            "$$\\text{RSS} = \\sum_i (y_i - \\hat{y}_i)^2$$",
            "square every residual and add them up. This 'residual sum of squares' measures how badly",
            "the line misses overall. The ", strong("least-squares line"), " is the one with the smallest RSS.",
            "Squaring stops positive and negative misses from cancelling out, and it punishes big",
            "misses much more than small ones. Turn on 'Show squared residuals' to see the squares",
            "whose total area you're trying to shrink.",
            symbols = list("\\(y_i\\)" = "actual value of point i", "\\(\\hat{y}_i\\)" = "its prediction")
          )
        )),
        layout_columns(
          col_widths = c(7, 5),
          card(
            full_screen = TRUE,
            card_header("The data and your line"),
            plotOutput(ns("fit_plot"), height = "470px"),
            axes_note("x, the predictor (the variable we use to predict).",
                      "y, the response (the variable we're trying to predict).",
                      "Orange: your line. Gray segments: residuals. Green dashed: the least-squares line, once revealed.")
          ),
          div(
            uiOutput(ns("fit_stats")),
            card(
              card_header("The RSS landscape"),
              plotOutput(ns("rss_map"), height = "300px"),
              axes_note("intercept of a candidate line.",
                        "slope of a candidate line.",
                        "Every spot on this map is a possible line; darker = larger RSS (worse fit).",
                        "The orange dot is your line. Least squares finds the lowest point of the valley (green).")
            )
          )
        )
      ),

      # ==== 2. Formulas & R² ===================================================
      nav_panel(
        "2 · Formulas & R²", value = "formula",
        explain(p(
          "You don't need to search for the best line by trial and error. A little algebra gives",
          "the answer directly. Below, each formula comes with the numbers from the current sample plugged in."
        )),
        layout_columns(
          col_widths = c(6, 6),
          card(
            card_header("The slope and intercept"),
            withMathJax(
              math_words(
                "$$b_1 = \\frac{\\sum_i (x_i - \\bar{x})(y_i - \\bar{y})}{\\sum_i (x_i - \\bar{x})^2}$$",
                "how much x and y move together, divided by how much x moves on its own.",
                "If y tends to be above its average whenever x is above its average, the top is",
                "positive and so is the slope.",
                symbols = list("\\(\\bar{x}, \\bar{y}\\)" = "the averages of x and y")
              ),
              math_words(
                "$$b_0 = \\bar{y} - b_1\\,\\bar{x}$$",
                "choose the intercept so that the line passes through the 'point of averages' (x̄, ȳ).",
                "The least-squares line always goes through the middle of the data."
              ),
              deeper(
                "Show where these come from",
                p("RSS is a bowl-shaped function of \\(b_0\\) and \\(b_1\\) (the landscape in tab 1).",
                  "At the bottom of a bowl the slope is flat in every direction, so set both derivatives to zero:",
                  "$$\\frac{\\partial \\text{RSS}}{\\partial b_0} = -2\\sum_i (y_i - b_0 - b_1 x_i) = 0,",
                  "\\qquad \\frac{\\partial \\text{RSS}}{\\partial b_1} = -2\\sum_i x_i (y_i - b_0 - b_1 x_i) = 0$$",
                  "Solving these two equations together gives the formulas above.")
              ),
              uiOutput(ns("formula_numbers"))
            )
          ),
          card(
            card_header("How good is the line? R²"),
            withMathJax(
              p("Start with how much y varies around its own average. The line 'explains' part",
                "of that variation; the residuals are the part it doesn't."),
              math_words(
                "$$\\underset{\\text{total}}{\\text{SST}} = \\underset{\\text{explained}}{\\text{SSR}} + \\underset{\\text{left over}}{\\text{SSE}}$$",
                "the total spread of y (squared distances from its average) splits exactly into",
                "the part captured by the line and the part left in the residuals. It's the same",
                "'kept + lost' idea as in PCA."
              ),
              math_words(
                "$$R^2 = \\frac{\\text{SSR}}{\\text{SST}} = 1 - \\frac{\\text{SSE}}{\\text{SST}}$$",
                "the fraction of y's variation that the line explains: 0 = none, 1 = all of it",
                "(every point exactly on the line)."
              )
            ),
            uiOutput(ns("r2_box"))
          )
        ),
        card(
          card_header("The three sums of squares, drawn"),
          plotOutput(ns("ss_plot"), height = "330px"),
          axes_note("x, the predictor.", "y, the response.",
                    "Left: each point's distance from the average of y (dotted). Middle: how far the line moves",
                    "each prediction away from that average. Right: the residuals. Square the segments in each",
                    "panel and add them up to get SST, SSR and SSE.")
        )
      ),

      # ==== 3. How sure are we? ================================================
      nav_panel(
        "3 · How sure are we?", value = "uncert",
        explain(withMathJax(
          p("Our data are just one random sample. A different sample would give a slightly different",
            "line. To see how much the line wobbles, the app redraws 200 new samples from the same",
            "'true' line and fits each one."),
          math_words(
            "$$\\text{SE}(b_1) = \\frac{\\hat{\\sigma}}{\\sqrt{\\sum_i (x_i - \\bar{x})^2}}$$",
            "the ", strong("standard error"), " is the typical size of the slope's wobble from sample to sample.",
            "It shrinks when there's less noise (smaller σ̂), more data, or more spread-out x values",
            "(a bigger bottom). Try the sliders.",
            symbols = list("\\(\\hat{\\sigma}\\)" = "estimated noise: the typical size of a residual")
          ),
          p("From the standard error we get a ", strong("95% confidence interval"),
            " (slope ± about 2 SEs), a range that captures the true slope in 95% of samples.",
            "We also get a ", strong("t-statistic"), " (slope ÷ SE: how many SEs the slope is from zero)",
            "and a ", strong("p-value"), " (how surprising a slope this far from zero would be if the",
            "true slope were 0).")
        )),
        layout_columns(
          col_widths = c(6, 6),
          card(full_screen = TRUE,
               card_header("This sample, with uncertainty bands"),
               plotOutput(ns("band_plot"), height = "380px"),
               axes_note("x, the predictor.", "y, the response.",
                         "Blue line: fitted to this sample. Green dashed: the true line. Bands: 95% ranges.")),
          card(full_screen = TRUE,
               card_header("200 other samples, 200 other lines"),
               plotOutput(ns("spaghetti"), height = "380px"),
               axes_note("x, the predictor.", "y, the response.",
                         "Each thin gray line is the least-squares line from one simulated sample.",
                         "They fan out most at the ends, where small changes in slope matter most."))
        ),
        layout_columns(
          col_widths = c(6, 6),
          card(card_header("The slope across samples"),
               plotOutput(ns("slope_hist"), height = "300px"),
               axes_note("estimated slope.", "how many of the 200 samples gave that slope.",
                         "The curve is what the standard-error formula predicts; the bars are what actually",
                         "happened in the simulation. Green line: the true slope.")),
          card(card_header("Inference for this sample"),
               uiOutput(ns("infer_table")))
        )
      ),

      # ==== 4. Assumptions =====================================================
      nav_panel(
        "4 · Checking assumptions", value = "assume",
        explain(
          p("The least-squares line can always be computed, but the standard errors, intervals and",
            "p-values from tab 3 are only trustworthy if some assumptions hold. A handy mnemonic is ",
            strong("LINE"), ":"),
          tags$ul(
            tags$li(strong("L"), "inearity: the average of y really does follow a straight line in x."),
            tags$li(strong("I"), "ndependence: one point's error tells you nothing about another's",
                    "(violated by, e.g., repeated measures on the same people; see the mixed-effects page)."),
            tags$li(strong("N"), "ormality: the errors follow a bell curve (matters most for small samples)."),
            tags$li(strong("E"), "qual variance: the scatter around the line is the same size everywhere.")
          ),
          p("The errors themselves are invisible, so we inspect the residuals instead. Pick a scenario in",
            "the sidebar and compare the plots with 'Everything fine'.")
        ),
        layout_columns(
          col_widths = c(4, 4, 4),
          card(card_header("Data and fitted line"),
               plotOutput(ns("a_fit"), height = "300px"),
               axes_note("x, the predictor.", "y, the response.")),
          card(card_header("Residuals vs fitted values"),
               plotOutput(ns("a_resid"), height = "300px"),
               axes_note("the fitted (predicted) value for each point.",
                         "its residual (actual minus predicted).",
                         "Want: a shapeless band around 0. The orange curve traces the average residual;",
                         "it should stay flat.")),
          card(card_header("Normal Q-Q plot"),
               plotOutput(ns("a_qq"), height = "300px"),
               axes_note("where each residual would sit if the errors were perfectly bell-shaped.",
                         "where it actually sits (sorted residuals, standardized).",
                         "Want: points hugging the diagonal line. Curves at the ends signal skew or heavy tails."))
        )
      ),

      # ==== 5. Outliers & influence ============================================
      nav_panel(
        "5 · Outliers & influence", value = "influence",
        explain(
          p("Not all unusual points are equally dangerous. Two ideas help:"),
          tags$ul(
            tags$li(strong("Leverage"), " measures how unusual a point's x value is. Points far from the",
                    "middle of the x range act like a long lever: they can tilt the line a lot."),
            tags$li(strong("Influence"), " is how much the fitted line actually changes because of the point.",
                    "A point is influential when it has high leverage ", em("and"), " sits far from the",
                    "line the other points suggest. A common measure is ", strong("Cook's distance"),
                    "; values above about 1 deserve a closer look.")
          ),
          p("Try adding a point far above the line in the middle of the x range, then one far to",
            "the right. Which moves the line more?")
        ),
        layout_columns(
          col_widths = c(8, 4),
          card(full_screen = TRUE,
               card_header("Click to add points"),
               plotOutput(ns("infl_plot"), height = "480px", click = ns("infl_click")),
               axes_note("x, the predictor.", "y, the response.",
                         "Gray dashed: the line fitted to the original 15 points only. Blue: the line with your",
                         "points (orange) included.")),
          card(card_header("Your added points"),
               uiOutput(ns("infl_summary")),
               tableOutput(ns("infl_table")))
        )
      ),

      nav_panel("Quiz", value = "quiz", quiz_panel(ns))
    )
  )
}

linreg_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    B0_TRUE <- 2

    # ---- simulated sample shared by tabs 1-3
    seed <- reactiveVal(11)
    observeEvent(input$resample, seed(sample.int(1e6, 1)))
    dat <- reactive(sim_linear(input$n, B0_TRUE, input$b1_true, input$sd, seed()))
    fit <- reactive(lm(y ~ x, data = dat()))
    ols <- reactive(unname(coef(fit())))

    revealed <- reactiveVal(FALSE)
    observeEvent(list(input$n, input$b1_true, input$sd, input$resample), revealed(FALSE),
                 ignoreInit = TRUE)
    observeEvent(input$reveal, {
      revealed(TRUE)
      updateSliderInput(session, "my_b0", value = round(ols()[1], 1))
      updateSliderInput(session, "my_b1", value = round(ols()[2] / 0.05) * 0.05)
    })

    rss_of <- function(d, b0, b1) sum((d$y - b0 - b1 * d$x)^2)

    # ---- 1. the best line ---------------------------------------------------
    output$fit_plot <- renderPlot({
      d <- dat()
      b0 <- input$my_b0; b1 <- input$my_b1
      d$yhat <- b0 + b1 * d$x
      d$r <- d$y - d$yhat
      xl <- c(0, 10)
      yl <- range(c(d$y, d$yhat, B0_TRUE + c(0, 10) * input$b1_true)) + c(-1, 1)
      g <- ggplot(d)
      if (isTRUE(input$show_sq)) {
        w <- abs(d$r) * diff(xl) / diff(yl)    # visual squares on an aspect-1 panel
        side <- ifelse(d$x + w > xl[2], -1, 1)
        g <- g + geom_rect(aes(xmin = x, xmax = x + side * w, ymin = pmin(y, yhat), ymax = pmax(y, yhat)),
                           fill = sv_col$accent, alpha = 0.12, colour = sv_col$accent, linewidth = 0.2)
      }
      if (isTRUE(input$show_resid)) {
        g <- g + geom_segment(aes(x, y, xend = x, yend = yhat), colour = sv_col$muted, linewidth = 0.5)
      }
      if (revealed()) {
        g <- g + geom_abline(intercept = ols()[1], slope = ols()[2], colour = sv_col$good,
                             linewidth = 1, linetype = "22")
      }
      g + geom_abline(intercept = b0, slope = b1, colour = sv_col$accent, linewidth = 1.2) +
        geom_point(aes(x, y), colour = sv_col$primary, size = 2.6, alpha = 0.85) +
        coord_cartesian(xlim = xl, ylim = yl, expand = FALSE) +
        labs(x = "x", y = "y") +
        theme_statviz() + theme(aspect.ratio = 1)
    })

    output$fit_stats <- renderUI({
      d <- dat()
      mine <- rss_of(d, input$my_b0, input$my_b1)
      best <- rss_of(d, ols()[1], ols()[2])
      layout_column_wrap(
        width = 1 / 2, fill = FALSE,
        sv_box("Your RSS", fmt(mine, 1), p(sprintf("%.1f× the best", mine / best)), tone = "orange"),
        sv_box(if (revealed()) "Least-squares RSS" else "Best possible",
               if (revealed()) fmt(best, 1) else "?",
               p(if (revealed()) sprintf("intercept %s, slope %s", fmt(ols()[1]), fmt(ols()[2]))
                 else "press 'Show the least-squares line'"),
               tone = "green")
      )
    })

    output$rss_map <- renderPlot({
      d <- dat()
      grid <- expand.grid(b0 = seq(-5, 12, length.out = 90), b1 = seq(-2, 3, length.out = 90))
      sy <- sum(d$y^2); sx <- sum(d$x); sx2 <- sum(d$x^2); sY <- sum(d$y); sxy <- sum(d$x * d$y)
      n <- nrow(d)
      grid$rss <- with(grid, sy + n * b0^2 + b1^2 * sx2 - 2 * b0 * sY - 2 * b1 * sxy + 2 * b0 * b1 * sx)
      g <- ggplot(grid, aes(b0, b1)) +
        geom_raster(aes(fill = log10(rss)), interpolate = TRUE) +
        geom_contour(aes(z = log10(rss)), colour = "#ffffff", alpha = 0.5, bins = 12, linewidth = 0.3) +
        scale_fill_gradientn(colours = seq_blues, name = "RSS",
                             labels = function(v) format(round(10^v), big.mark = ",")) +
        annotate("point", x = input$my_b0, y = input$my_b1, colour = sv_col$surface, size = 5.5) +
        annotate("point", x = input$my_b0, y = input$my_b1, colour = sv_col$accent, size = 4)
      if (revealed()) {
        g <- g + annotate("point", x = ols()[1], y = ols()[2], colour = sv_col$surface, size = 5.5) +
          annotate("point", x = ols()[1], y = ols()[2], colour = sv_col$good, size = 4, shape = 18)
      }
      g + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
        labs(x = "Intercept", y = "Slope") +
        theme_statviz() + theme(legend.position = "right")
    })

    # ---- 2. formulas & R² ---------------------------------------------------
    output$formula_numbers <- renderUI({
      d <- dat()
      xb <- mean(d$x); yb <- mean(d$y)
      sxy <- sum((d$x - xb) * (d$y - yb)); sxx <- sum((d$x - xb)^2)
      b1 <- sxy / sxx; b0 <- yb - b1 * xb
      withMathJax(div(
        class = "sv-math",
        p(strong("With this sample's numbers:")),
        sprintf("$$b_1 = \\frac{%s}{%s} = %s \\qquad b_0 = %s - %s \\times %s = %s$$",
                fmt(sxy, 1), fmt(sxx, 1), fmt(b1), fmt(yb), fmt(b1), fmt(xb), fmt(b0)),
        p(class = "mb-0", "So the fitted line is ", strong(sprintf("ŷ = %s + %s x", fmt(b0), fmt(b1))), ".",
          "In words: each 1-unit increase in x goes with a ", fmt(abs(b1)),
          if (b1 >= 0) " increase" else " decrease", " in predicted y.",
          sprintf("(The true slope used to make the data is %s.)", fmt(input$b1_true)))
      ))
    })

    ss <- reactive({
      d <- dat(); yhat <- fitted(fit()); yb <- mean(d$y)
      c(SST = sum((d$y - yb)^2), SSR = sum((yhat - yb)^2), SSE = sum((d$y - yhat)^2))
    })

    output$r2_box <- renderUI({
      s <- ss()
      layout_column_wrap(
        width = 1 / 2, fill = FALSE,
        sv_box("R²", sprintf("%.0f%%", 100 * s["SSR"] / s["SST"]),
               p(sprintf("SSR %s ÷ SST %s", fmt(s["SSR"], 1), fmt(s["SST"], 1))), tone = "blue"),
        sv_box("Left unexplained", sprintf("%.0f%%", 100 * s["SSE"] / s["SST"]),
               p(sprintf("SSE %s", fmt(s["SSE"], 1))), tone = "gray")
      )
    })

    output$ss_plot <- renderPlot({
      d <- dat(); yhat <- fitted(fit()); yb <- mean(d$y); s <- ss()
      panels <- c(sprintf("Total: y − ȳ   (SST = %s)", fmt(s["SST"], 1)),
                  sprintf("Explained: ŷ − ȳ   (SSR = %s)", fmt(s["SSR"], 1)),
                  sprintf("Left over: y − ŷ   (SSE = %s)", fmt(s["SSE"], 1)))
      seg <- rbind(
        data.frame(panel = panels[1], x = d$x, y = d$y, yend = yb),
        data.frame(panel = panels[2], x = d$x, y = yhat, yend = yb),
        data.frame(panel = panels[3], x = d$x, y = d$y, yend = yhat)
      )
      seg$panel <- factor(seg$panel, levels = panels)
      pts <- data.frame(panel = factor(panels, levels = panels))
      ggplot(seg) +
        geom_hline(yintercept = yb, colour = sv_col$ink2, linetype = "dotted") +
        geom_abline(data = pts[2:3, , drop = FALSE], aes(intercept = coef(fit())[1], slope = coef(fit())[2]),
                    colour = sv_col$primary, linewidth = 0.9) +
        geom_segment(aes(x, y, xend = x, yend = yend), colour = sv_col$accent, linewidth = 0.6) +
        geom_point(data = do.call(rbind, lapply(panels, function(p) data.frame(panel = p, x = d$x, y = d$y))) |>
                     transform(panel = factor(panel, levels = panels)),
                   aes(x, y), colour = sv_col$primary, alpha = 0.6, size = 1.8) +
        facet_wrap(~panel) +
        labs(x = "x", y = "y") +
        theme_statviz() + theme(strip.text = element_text(face = "bold", colour = sv_col$ink))
    })

    # ---- 3. uncertainty -----------------------------------------------------
    reps <- reactive({
      n <- input$n; b1 <- input$b1_true; sd <- input$sd; s0 <- seed()
      t(vapply(seq_len(200), function(i) {
        d <- sim_linear(n, B0_TRUE, b1, sd, s0 * 7 + i)
        sxx <- sum((d$x - mean(d$x))^2)
        slope <- sum((d$x - mean(d$x)) * (d$y - mean(d$y))) / sxx
        icpt <- mean(d$y) - slope * mean(d$x)
        sig <- sqrt(sum((d$y - icpt - slope * d$x)^2) / (n - 2))
        se <- sig / sqrt(sxx)
        tq <- qt(0.975, n - 2)
        c(b0 = icpt, b1 = slope, covered = abs(slope - b1) <= tq * se)
      }, numeric(3)))
    })

    output$band_plot <- renderPlot({
      d <- dat(); f <- fit()
      grid <- data.frame(x = seq(0, 10, length.out = 100))
      ci <- predict(f, grid, interval = "confidence")
      pi <- predict(f, grid, interval = "prediction")
      g <- ggplot()
      if (isTRUE(input$show_pi)) {
        g <- g + geom_ribbon(aes(x = grid$x, ymin = pi[, "lwr"], ymax = pi[, "upr"]),
                             fill = sv_col$primary, alpha = 0.12)
      }
      if (isTRUE(input$show_ci)) {
        g <- g + geom_ribbon(aes(x = grid$x, ymin = ci[, "lwr"], ymax = ci[, "upr"]),
                             fill = sv_col$primary, alpha = 0.3)
      }
      g + geom_abline(intercept = B0_TRUE, slope = input$b1_true, colour = sv_col$good,
                      linetype = "22", linewidth = 0.9) +
        geom_line(aes(grid$x, ci[, "fit"]), colour = sv_col$primary, linewidth = 1.1) +
        geom_point(data = d, aes(x, y), colour = sv_col$primary, alpha = 0.7, size = 2.2) +
        labs(x = "x", y = "y") + theme_statviz()
    })

    output$spaghetti <- renderPlot({
      r <- as.data.frame(reps())
      d <- dat()
      ggplot() +
        geom_abline(data = r, aes(intercept = b0, slope = b1), colour = sv_col$muted, alpha = 0.18) +
        geom_abline(intercept = B0_TRUE, slope = input$b1_true, colour = sv_col$good,
                    linetype = "22", linewidth = 1) +
        geom_point(data = d, aes(x, y), alpha = 0) +
        coord_cartesian(xlim = c(0, 10)) +
        labs(x = "x", y = "y") + theme_statviz()
    })

    output$slope_hist <- renderPlot({
      r <- as.data.frame(reps())
      n <- input$n
      se_theory <- input$sd / sqrt((n - 1) * 100 / 12)    # Var(Uniform(0, 10)) = 100/12
      bw <- diff(range(r$b1)) / 25
      curve <- data.frame(x = seq(min(r$b1), max(r$b1), length.out = 200))
      curve$y <- dnorm(curve$x, input$b1_true, se_theory) * nrow(r) * bw
      ggplot(r, aes(b1)) +
        geom_histogram(binwidth = bw, fill = "#9ec5f4", colour = sv_col$surface) +
        geom_line(data = curve, aes(x, y), colour = sv_col$primary, linewidth = 1) +
        geom_vline(xintercept = input$b1_true, colour = sv_col$good, linewidth = 1, linetype = "22") +
        labs(x = "Estimated slope", y = "Number of samples") + theme_statviz()
    })

    output$infer_table <- renderUI({
      f <- summary(fit())$coefficients
      ci <- confint(fit())
      cover <- mean(reps()[, "covered"])
      b <- f["x", ]
      tagList(
        tags$table(
          class = "table table-sm",
          tags$tbody(
            tags$tr(tags$td("Estimated slope"), tags$td(strong(fmt(b["Estimate"], 3)))),
            tags$tr(tags$td("Standard error"), tags$td(fmt(b["Std. Error"], 3))),
            tags$tr(tags$td("95% confidence interval"),
                    tags$td(sprintf("%s to %s", fmt(ci["x", 1], 3), fmt(ci["x", 2], 3)))),
            tags$tr(tags$td("t = slope ÷ SE"), tags$td(fmt(b["t value"]))),
            tags$tr(tags$td("p-value"), tags$td(fmt_p(b["Pr(>|t|)"])))
          )
        ),
        p(class = "small",
          sprintf("In words: our best guess for the slope is %s, and we're fairly confident the true slope is between %s and %s.",
                  fmt(b["Estimate"]), fmt(ci["x", 1]), fmt(ci["x", 2])),
          if (b["Pr(>|t|)"] < 0.05) "The interval excludes 0, so the data give clear evidence of a relationship."
          else "The interval includes 0, so we can't rule out 'no relationship at all'."),
        div(class = "alert alert-secondary py-2 small mb-0",
            icon("circle-check"),
            sprintf(" Of the 200 simulated samples, %.0f%% of their 95%% intervals contained the true slope (%s).",
                    100 * cover, fmt(input$b1_true)),
            "That's what '95% confidence' means: the method works 95% of the time.")
      )
    })

    # ---- 4. assumptions -----------------------------------------------------
    seed_a <- reactiveVal(5)
    observeEvent(input$resample_a, seed_a(sample.int(1e6, 1)))
    dat_a <- reactive(sim_linear(input$n_a, B0_TRUE, 0.8, 1.5, seed_a(), input$scenario))
    fit_a <- reactive(lm(y ~ x, data = dat_a()))

    output$a_fit <- renderPlot({
      ggplot(dat_a(), aes(x, y)) +
        geom_point(colour = sv_col$primary, alpha = 0.7, size = 2) +
        geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = sv_col$ink, linewidth = 0.9) +
        theme_statviz()
    })

    output$a_resid <- renderPlot({
      f <- fit_a()
      df <- data.frame(fitted = fitted(f), resid = resid(f))
      ggplot(df, aes(fitted, resid)) +
        geom_hline(yintercept = 0, colour = sv_col$ink2) +
        geom_point(colour = sv_col$primary, alpha = 0.7, size = 2) +
        geom_smooth(method = "loess", formula = y ~ x, se = FALSE, colour = sv_col$accent, linewidth = 1) +
        labs(x = "Fitted value", y = "Residual") + theme_statviz()
    })

    output$a_qq <- renderPlot({
      df <- data.frame(r = rstandard(fit_a()))
      ggplot(df, aes(sample = r)) +
        stat_qq_line(colour = sv_col$ink2) +
        stat_qq(colour = sv_col$primary, alpha = 0.7, size = 2) +
        labs(x = "Theoretical (bell-curve) quantiles", y = "Standardized residuals") + theme_statviz()
    })

    output$assume_advice <- renderUI({
      txt <- switch(
        input$scenario,
        ok = "All four assumptions hold. Residuals form a flat, even band, and the Q-Q points follow the line. Use this as your reference picture.",
        curve = "Linearity is broken: the residual plot shows a U-shape (the orange curve bends). The straight line systematically over- and under-predicts. Fix: add a curved term (e.g. x²) or transform a variable.",
        hetero = "Equal variance is broken: the residuals fan out like a funnel. The line itself is fine, but the standard errors and p-values are off. Fix: transform y (e.g. log) or use robust standard errors.",
        skew = "Normality is broken: the Q-Q plot curves away from the line, because a few residuals are much larger than a bell curve allows. With large samples this matters little; with small samples, intervals can mislead.",
        outlier = "One point sits far from the rest: it stands out in all three plots and drags the line toward itself. Check whether it's a data error before deciding what to do. See tab 5."
      )
      div(class = "alert alert-info small mt-3", strong("What to look for: "), txt)
    })

    # ---- 5. influence -------------------------------------------------------
    base <- sim_linear(15, B0_TRUE, 0.8, 1.2, 3)
    added <- reactiveVal(data.frame(x = numeric(0), y = numeric(0)))
    observeEvent(input$infl_click, {
      cl <- input$infl_click
      added(rbind(added(), data.frame(x = cl$x, y = cl$y)))
    })
    observeEvent(input$undo, { a <- added(); if (nrow(a)) added(a[-nrow(a), , drop = FALSE]) })
    observeEvent(input$reset, added(data.frame(x = numeric(0), y = numeric(0))))

    infl_fit <- reactive(lm(y ~ x, data = rbind(base, added())))

    output$infl_plot <- renderPlot({
      a <- added()
      f0 <- coef(lm(y ~ x, data = base)); f1 <- coef(infl_fit())
      g <- ggplot() +
        geom_abline(intercept = f0[1], slope = f0[2], colour = sv_col$muted, linetype = "22", linewidth = 1) +
        geom_abline(intercept = f1[1], slope = f1[2], colour = sv_col$primary, linewidth = 1.2) +
        geom_point(data = base, aes(x, y), colour = sv_col$primary, size = 2.8, alpha = 0.8)
      if (nrow(a)) {
        a$lab <- seq_len(nrow(a))
        g <- g + geom_point(data = a, aes(x, y), colour = sv_col$accent, size = 4) +
          geom_text(data = a, aes(x, y, label = lab), colour = "#ffffff", size = 2.8, fontface = "bold")
      }
      g + coord_cartesian(xlim = c(-2, 20), ylim = c(-10, 25), expand = FALSE) +
        labs(x = "x", y = "y") + theme_statviz()
    })

    output$infl_summary <- renderUI({
      f0 <- coef(lm(y ~ x, data = base)); f1 <- coef(infl_fit())
      tagList(
        p(sprintf("Slope without your points: %s", fmt(f0[2]))),
        p("Slope with them: ", strong(fmt(f1[2])),
          sprintf(" (%s%s)", if (f1[2] >= f0[2]) "+" else "", fmt(f1[2] - f0[2]))),
        if (!nrow(added())) p(class = "text-muted small", "No points added yet.")
      )
    })

    output$infl_table <- renderTable({
      a <- added()
      req(nrow(a) > 0)
      f <- infl_fit()
      idx <- nrow(base) + seq_len(nrow(a))
      data.frame(`#` = seq_len(nrow(a)), x = a$x, y = a$y,
                 Leverage = hatvalues(f)[idx], `Cook's D` = cooks.distance(f)[idx],
                 check.names = FALSE)
    }, digits = 2, striped = TRUE, spacing = "s")

    # ---- quiz ---------------------------------------------------------------
    quiz_server(input, output, session, list(
      list(q = "A residual is...",
           choices = c("the slope of the line", "the vertical gap between a point and the line (actual − predicted)",
                       "the perpendicular distance from a point to the line", "the average of y"),
           answer = 2,
           why = "Regression measures misses vertically, because it's predicting y. (PCA uses perpendicular distances.)"),
      list(q = "Least squares chooses the line that...",
           choices = c("passes through the most points", "makes the sum of the residuals zero and nothing else",
                       "minimizes the sum of squared residuals", "has the steepest slope"),
           answer = 3,
           why = "It minimizes RSS, the total area of the squares in tab 1."),
      list(q = "R² = 0.8 means...",
           choices = c("80% of the points are on the line", "the slope is 0.8",
                       "the line explains 80% of the variation in y", "the prediction is right 80% of the time"),
           answer = 3,
           why = "R² is SSR / SST: the share of y's total variation captured by the line."),
      list(q = "Which change makes the slope's standard error SMALLER?",
           choices = c("More noise", "Fewer data points", "x values bunched close together",
                       "x values spread out over a wider range"),
           answer = 4,
           why = "SE = σ̂ / √Σ(x − x̄)²: spreading x out makes the bottom bigger, so the slope is pinned down more precisely."),
      list(q = "The residuals-vs-fitted plot shows a clear U-shape. Which assumption is in trouble?",
           choices = c("Linearity", "Independence", "Normality", "Equal variance"),
           answer = 1,
           why = "A curve in the residuals means a straight line isn't capturing the real shape of the relationship."),
      list(q = "Which added point would most change the fitted line?",
           choices = c("One near the middle of the x range, close to the line",
                       "One near the middle of the x range, far above the line",
                       "One far beyond the other x values, far from the line",
                       "One far beyond the other x values, right on the line"),
           answer = 3,
           why = "High leverage (unusual x) plus a big residual = high influence.")
    ))
  })
}
