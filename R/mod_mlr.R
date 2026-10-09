# Multiple linear regression module --------------------------------------------
#
# Tabs:
#   1. A plane, not a line     - 3-D view of two predictors
#   2. Holding others constant - confounding; simple vs multiple slope; added-variable plot
#   3. Multicollinearity       - correlated predictors inflate uncertainty
#   4. Real data               - coefficients, intervals, plain-language interpretation
#   5. Quiz

mlr_ui <- function(id) {
  ns <- NS(id)

  layout_sidebar(
    sidebar = sidebar(
      width = 300,
      show_on(ns, c("plane", "adjust", "collin"),
        h6("Simulated data"),
        p(class = "small text-muted",
          "y is built from two predictors, x1 and x2, plus noise. You set the true effects."),
        sliderInput(ns("n"), "Number of points", 20, 300, 80, step = 10),
        sliderInput(ns("b1"), "True effect of x1", -2, 2, 0.5, step = 0.1),
        sliderInput(ns("b2"), "True effect of x2", -2, 2, 1.5, step = 0.1),
        sliderInput(ns("rho"), "Correlation between x1 and x2", -0.95, 0.95, 0.7, step = 0.05),
        sliderInput(ns("sd"), "Noise (SD)", 0.5, 5, 1.5, step = 0.1),
        actionButton(ns("resample"), "Draw a new sample", icon = icon("shuffle"),
                     class = "btn-outline-secondary btn-sm w-100")
      ),
      show_on(ns, "plane",
        hr(),
        checkboxInput(ns("show_plane"), "Show fitted plane", TRUE),
        checkboxInput(ns("show_resid3"), "Show residuals", TRUE)
      ),
      show_on(ns, "real",
        selectInput(ns("dataset"), "Dataset",
                    setNames(names(mlr_datasets), vapply(mlr_datasets, `[[`, "", "name"))),
        uiOutput(ns("blurb")),
        selectInput(ns("yvar"), "Response (what we predict)", NULL),
        selectizeInput(ns("xvars"), "Predictors", NULL, multiple = TRUE,
                       options = list(plugins = list("remove_button"))),
        checkboxInput(ns("standardize"), "Standardize predictors (compare effect sizes)", FALSE)
      ),
      show_on(ns, "quiz", p("Test yourself on the ideas from the other tabs."))
    ),

    navset_underline(
      id = ns("tab"),

      # ==== 1. Plane ===========================================================
      nav_panel(
        "1 · A plane, not a line", value = "plane",
        explain(withMathJax(
          p("With two predictors, each observation has three numbers (x1, x2, y), so the data form a",
            "cloud in 3-D. Instead of a line, regression fits a flat ", strong("plane"), " through it."),
          math_words(
            "$$\\hat{y} = b_0 + b_1\\,x_1 + b_2\\,x_2$$",
            "start at the intercept, add b1 for every unit of x1 and b2 for every unit of x2.",
            "The idea is the same as before: choose the b's that make the sum of squared residuals",
            "(vertical gaps between the points and the plane) as small as possible.",
            symbols = list("\\(b_1\\)" = "tilt of the plane in the x1 direction",
                           "\\(b_2\\)" = "tilt in the x2 direction")
          ),
          p("Drag the plot to rotate it. Look along the x2 axis and the plane's edge shows you b1;",
            "look along x1 and you see b2. With more than two predictors it's a 'hyperplane' we can't",
            "draw, but the math is identical.")
        )),
        layout_columns(
          col_widths = c(8, 4),
          card(full_screen = TRUE, card_header("Data and fitted plane"),
               plotlyOutput(ns("plane3d"), height = "540px"),
               card_footer(class = "small sv-axes",
                           "Axes: x1 and x2 (the predictors) along the floor, y (the response) going up.",
                           "Gray lines: residuals.")),
          div(uiOutput(ns("plane_stats")))
        )
      ),

      # ==== 2. Adjusting =======================================================
      nav_panel(
        "2 · Holding others constant", value = "adjust",
        explain(withMathJax(
          p("Each coefficient in a multiple regression answers: ",
            em("'how much does y change when this predictor goes up by 1, ",
               strong("holding the other predictors fixed"), "?'"),
            "That can be very different from the slope you get using one predictor alone."),
          p(strong("Why?"), "When x1 and x2 are correlated, x1 'carries' some of x2's effect.",
            "A simple regression of y on x1 alone gives x1 credit for both (orange line below).",
            "Multiple regression separates them. Try: set the true effect of x1 to 0, keep x2's",
            "effect large and the correlation high. The simple slope still looks strong,",
            "but the multiple-regression coefficient is close to 0. This is ", strong("confounding"), "."),
          math_words(
            "$$\\text{simple slope} \\;\\approx\\; b_1 + b_2 \\times (\\text{slope of } x_2 \\text{ on } x_1)$$",
            "leaving x2 out doesn't make its effect disappear; it gets absorbed into x1's slope,",
            "in proportion to how strongly x2 tracks x1. This is called ", strong("omitted-variable bias"), "."
          )
        )),
        layout_columns(
          col_widths = c(6, 6),
          card(full_screen = TRUE, card_header("y against x1, shaded by x2"),
               plotOutput(ns("adjust_plot"), height = "400px"),
               axes_note("x1.", "y.",
                         "Point shade = value of x2 (darker = higher). Orange: the simple regression of y on x1.",
                         "Blue: the multiple-regression slope for x1, drawn at low, middle and high x2. They're parallel:",
                         "within each slice of x2, y rises by b1 per unit of x1.")),
          card(full_screen = TRUE, card_header("Added-variable plot: x1's effect after removing x2"),
               plotOutput(ns("avp_plot"), height = "400px"),
               axes_note("the part of x1 that x2 can't predict (residuals of x1 regressed on x2).",
                         "the part of y that x2 can't predict (residuals of y regressed on x2).",
                         "Remove x2's influence from both, and the slope of what's left is exactly b1",
                         "from the multiple regression. This is what 'holding x2 constant' means."))
        ),
        uiOutput(ns("adjust_table"))
      ),

      # ==== 3. Collinearity ====================================================
      nav_panel(
        "3 · Multicollinearity", value = "collin",
        explain(withMathJax(
          p("If two predictors are strongly correlated, the data contain little information about",
            "what happens when one changes and the other doesn't. The model can still predict y well,",
            "but it struggles to say how to split the credit between them. The result:",
            "coefficients that wobble wildly from sample to sample."),
          math_words(
            "$$\\text{VIF} = \\frac{1}{1 - r^2}$$",
            "the ", strong("variance inflation factor"), " says how many times larger a coefficient's",
            "variance is because of the correlation, compared with uncorrelated predictors.",
            "VIF above about 5-10 is a common warning sign.",
            symbols = list("\\(r\\)" = "correlation between x1 and x2 (with more predictors: how well the others predict this one)")
          ),
          p("Slide the correlation toward ±0.95 and watch the cloud of estimates stretch out.")
        )),
        layout_columns(
          col_widths = c(7, 5),
          card(full_screen = TRUE, card_header("Estimates of b1 and b2 from 200 samples"),
               plotOutput(ns("collin_plot"), height = "430px"),
               axes_note("estimated effect of x1 (b1).", "estimated effect of x2 (b2).",
                         "Each dot is one simulated sample's pair of estimates. Green cross: the true values.",
                         "With correlated predictors the cloud tilts: when one estimate comes out too high,",
                         "the other tends to come out too low.")),
          div(uiOutput(ns("collin_stats")))
        )
      ),

      # ==== 4. Real data =======================================================
      nav_panel(
        "4 · Real data", value = "real",
        explain(withMathJax(
          p("Pick a dataset, a response and some predictors. The table translates each coefficient into",
            "a sentence. Check how the coefficients change as you add or remove predictors:",
            "that's 'holding others constant' at work."),
          deeper(
            "Show the matrix formula computers use",
            math_words(
              "$$\\mathbf{b} = (X^\\top X)^{-1} X^\\top \\mathbf{y}$$",
              "stack all the predictors as columns of a table X (plus a column of 1s for the intercept).",
              "This one formula solves 'minimize the sum of squared residuals' for any number of predictors at once.",
              "\\(X^\\top X\\) measures how the predictors overlap with each other; \\(X^\\top \\mathbf{y}\\) how",
              "each overlaps with y; the inverse 'divides out' the overlap between predictors.",
              symbols = list("\\(X^{-1}\\)" = "matrix inverse (the matrix version of dividing)")
            )
          )
        )),
        layout_columns(
          col_widths = c(6, 6),
          card(card_header("Coefficients with 95% confidence intervals"),
               plotOutput(ns("coef_plot"), height = "320px"),
               axes_note("the coefficient: change in the response per 1-unit increase in the predictor",
                         "the predictors.",
                         "Bars show 95% confidence intervals. If a bar crosses the dashed zero line, the data",
                         "are consistent with no effect (holding the others fixed). Predictors measured in different",
                         "units can't be compared by size; tick 'Standardize predictors' to put them on one scale.")),
          card(card_header("Model fit"),
               uiOutput(ns("fit_boxes")),
               plotOutput(ns("real_resid"), height = "220px"),
               axes_note("fitted (predicted) value.", "residual.",
                         "Want: a shapeless band around 0 (see 'Checking assumptions' on the simple regression page)."))
        ),
        card(card_header("What each coefficient means"), uiOutput(ns("coef_sentences")))
      ),

      nav_panel("Quiz", value = "quiz", quiz_panel(ns))
    )
  )
}

mlr_server <- function(id) {
  moduleServer(id, function(input, output, session) {

    seed <- reactiveVal(21)
    observeEvent(input$resample, seed(sample.int(1e6, 1)))
    dat <- reactive(sim_mlr(input$n, input$b1, input$b2, input$rho, input$sd, seed()))
    fit <- reactive(lm(y ~ x1 + x2, data = dat()))

    # ---- 1. plane -----------------------------------------------------------
    output$plane3d <- renderPlotly({
      d <- dat(); f <- fit(); b <- coef(f)
      p <- plot_ly() |>
        add_trace(data = d, x = ~x1, y = ~x2, z = ~y, type = "scatter3d", mode = "markers",
                  marker = list(size = 3.5, color = sv_col$primary), name = "data",
                  text = sprintf("x1 %.1f<br>x2 %.1f<br>y %.1f", d$x1, d$x2, d$y), hoverinfo = "text")
      if (isTRUE(input$show_resid3)) {
        yh <- fitted(f)
        xs <- as.vector(rbind(d$x1, d$x1, NA)); ys <- as.vector(rbind(d$x2, d$x2, NA))
        zs <- as.vector(rbind(d$y, yh, NA))
        p <- p |> add_trace(x = xs, y = ys, z = zs, type = "scatter3d", mode = "lines",
                            line = list(color = sv_col$muted, width = 2), hoverinfo = "skip", name = "residuals")
      }
      if (isTRUE(input$show_plane)) {
        g1 <- seq(min(d$x1), max(d$x1), length.out = 15); g2 <- seq(min(d$x2), max(d$x2), length.out = 15)
        z <- outer(g1, g2, function(a, c) b[1] + b[2] * a + b[3] * c)
        p <- p |> add_surface(x = g1, y = g2, z = t(z), opacity = 0.45, showscale = FALSE,
                              colorscale = list(c(0, 1), c(sv_col$good, sv_col$good)), hoverinfo = "skip")
      }
      p |> layout(scene = list(xaxis = list(title = "x1"), yaxis = list(title = "x2"), zaxis = list(title = "y")),
                  showlegend = FALSE, paper_bgcolor = sv_col$surface, margin = list(l = 0, r = 0, t = 0, b = 0)) |>
        config(displaylogo = FALSE)
    })

    output$plane_stats <- renderUI({
      b <- coef(fit())
      tagList(
        sv_box("Fitted effect of x1 (b1)", fmt(b["x1"]), p(sprintf("true value %s", fmt(input$b1))), tone = "blue"),
        sv_box("Fitted effect of x2 (b2)", fmt(b["x2"]), p(sprintf("true value %s", fmt(input$b2))), tone = "blue"),
        sv_box("R²", sprintf("%.0f%%", 100 * summary(fit())$r.squared),
               p("share of y's variation the plane explains"), tone = "green")
      )
    })

    # ---- 2. adjusting -------------------------------------------------------
    output$adjust_plot <- renderPlot({
      d <- dat(); b <- coef(fit()); s <- coef(lm(y ~ x1, data = d))
      q <- quantile(d$x2, c(0.1, 0.5, 0.9))
      slices <- data.frame(lvl = c("low x2", "middle x2", "high x2"), x2 = q,
                           a = b[1] + b[3] * q, col = c(seq_blues[2], seq_blues[4], seq_blues[7]))
      ggplot(d, aes(x1, y)) +
        geom_point(aes(colour = x2), size = 2.3, alpha = 0.9) +
        scale_colour_gradientn(colours = seq_blues[2:7], name = "x2") +
        geom_abline(data = slices, aes(intercept = a, slope = b[2]), colour = slices$col, linewidth = 1) +
        geom_abline(intercept = s[1], slope = s[2], colour = sv_col$accent, linewidth = 1.3) +
        annotate("label", x = max(d$x1), y = s[1] + s[2] * max(d$x1), label = sprintf("simple: %s", fmt(s[2])),
                 hjust = 1, vjust = -0.4, colour = sv_col$accent, fill = sv_col$surface, label.size = 0, fontface = "bold") +
        annotate("label", x = min(d$x1), y = slices$a[3] + b[2] * min(d$x1), label = sprintf("multiple: %s", fmt(b[2])),
                 hjust = 0, vjust = -0.4, colour = seq_blues[7], fill = sv_col$surface, label.size = 0, fontface = "bold") +
        theme_statviz() + theme(legend.position = "right")
    })

    output$avp_plot <- renderPlot({
      d <- dat()
      rx <- resid(lm(x1 ~ x2, data = d)); ry <- resid(lm(y ~ x2, data = d))
      ggplot(data.frame(rx, ry), aes(rx, ry)) +
        geom_hline(yintercept = 0, colour = sv_col$grid) + geom_vline(xintercept = 0, colour = sv_col$grid) +
        geom_point(colour = sv_col$primary, alpha = 0.7, size = 2.2) +
        geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = seq_blues[7], linewidth = 1.1) +
        labs(x = "x1, after removing x2", y = "y, after removing x2",
             subtitle = sprintf("slope = %s (same as b1)", fmt(coef(lm(ry ~ rx))[2]))) +
        theme_statviz()
    })

    output$adjust_table <- renderUI({
      d <- dat(); b <- coef(fit()); s <- coef(lm(y ~ x1, data = d))[2]
      g <- coef(lm(x2 ~ x1, data = d))[2]
      card(
        card_header("Comparing the two answers"),
        tags$table(
          class = "table table-sm mb-2",
          tags$thead(tags$tr(tags$th("Question"), tags$th("Model"), tags$th("Answer"))),
          tags$tbody(
            tags$tr(tags$td("How does y change with x1, ignoring x2?"), tags$td("y ~ x1"), tags$td(strong(fmt(s)))),
            tags$tr(tags$td("How does y change with x1, holding x2 fixed?"), tags$td("y ~ x1 + x2"),
                    tags$td(strong(fmt(b["x1"])))),
            tags$tr(tags$td("True effect of x1 (what we used to make the data)"), tags$td(""),
                    tags$td(fmt(input$b1)))
          )
        ),
        p(class = "small mb-0",
          sprintf("Check the bias formula: %s + %s × %s = %s, which matches the simple slope (%s).",
                  fmt(b["x1"]), fmt(b["x2"]), fmt(g), fmt(b["x1"] + b["x2"] * g), fmt(s)))
      )
    })

    # ---- 3. collinearity ----------------------------------------------------
    reps <- reactive({
      s0 <- seed()
      t(vapply(seq_len(200), function(i) {
        d <- sim_mlr(input$n, input$b1, input$b2, input$rho, input$sd, s0 * 5 + i)
        coef(lm(y ~ x1 + x2, data = d))[2:3]
      }, numeric(2)))
    })

    output$collin_plot <- renderPlot({
      r <- as.data.frame(reps()); names(r) <- c("b1", "b2")
      lim1 <- input$b1 + c(-1, 1) * max(abs(r$b1 - input$b1), abs(r$b2 - input$b2)) * 1.1
      lim2 <- input$b2 + c(-1, 1) * max(abs(r$b1 - input$b1), abs(r$b2 - input$b2)) * 1.1
      ggplot(r, aes(b1, b2)) +
        geom_point(colour = sv_col$primary, alpha = 0.5, size = 2) +
        annotate("point", x = input$b1, y = input$b2, colour = sv_col$good, shape = 3, size = 7, stroke = 2) +
        coord_equal(xlim = lim1, ylim = lim2) +
        labs(x = "Estimated b1", y = "Estimated b2") + theme_statviz()
    })

    output$collin_stats <- renderUI({
      d <- dat(); r <- reps()
      r_obs <- cor(d$x1, d$x2)
      se <- summary(fit())$coefficients[2:3, 2]
      tagList(
        sv_box("VIF (this sample)", fmt(1 / (1 - r_obs^2), 1),
               p(sprintf("correlation %s", fmt(r_obs))), tone = if (1 / (1 - r_obs^2) > 5) "orange" else "blue"),
        sv_box("Spread of b1 across samples", fmt(sd(r[, 1]), 3),
               p(sprintf("SE reported by this sample: %s", fmt(se[1], 3))), tone = "gray"),
        sv_box("R² (this sample)", sprintf("%.0f%%", 100 * summary(fit())$r.squared),
               p("prediction stays good even when the coefficients are shaky"), tone = "green")
      )
    })

    # ---- 4. real data -------------------------------------------------------
    ds <- reactive(mlr_datasets[[input$dataset]])
    output$blurb <- renderUI(p(class = "small text-muted", ds()$blurb))
    observeEvent(ds(), {
      nm <- names(ds()$data)
      updateSelectInput(session, "yvar", choices = nm, selected = ds()$y)
      updateSelectizeInput(session, "xvars", choices = setdiff(nm, ds()$y), selected = ds()$x)
    })
    observeEvent(input$yvar, {
      nm <- setdiff(names(ds()$data), input$yvar)
      updateSelectizeInput(session, "xvars", choices = nm, selected = intersect(input$xvars, nm))
    }, ignoreInit = TRUE)

    real_fit <- reactive({
      d <- ds()$data
      req(input$yvar %in% names(d))
      xv <- intersect(input$xvars, setdiff(names(d), input$yvar))
      validate(need(length(xv) >= 1, "Pick at least one predictor."))
      dd <- d[, c(input$yvar, xv), drop = FALSE]
      if (isTRUE(input$standardize)) dd[xv] <- lapply(dd[xv], function(v) as.numeric(scale(v)))
      lm(reformulate(sprintf("`%s`", xv), response = sprintf("`%s`", input$yvar)), data = dd)
    })

    coef_df <- reactive({
      f <- real_fit()
      s <- summary(f)$coefficients[-1, , drop = FALSE]
      ci <- confint(f)[-1, , drop = FALSE]
      data.frame(term = gsub("`", "", rownames(s)), est = s[, 1], se = s[, 2], p = s[, 4],
                 lo = ci[, 1], hi = ci[, 2])
    })

    output$coef_plot <- renderPlot({
      df <- coef_df()
      df$term <- factor(df$term, levels = rev(df$term))
      ggplot(df, aes(est, term)) +
        geom_vline(xintercept = 0, colour = sv_col$ink2, linetype = "22") +
        geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.2, colour = sv_col$primary, linewidth = 0.9) +
        geom_point(colour = sv_col$primary, size = 3.5) +
        labs(x = if (isTRUE(input$standardize)) "Change in response per 1 SD of predictor" else "Coefficient",
             y = NULL) +
        theme_statviz()
    })

    output$fit_boxes <- renderUI({
      s <- summary(real_fit())
      layout_column_wrap(
        width = 1 / 2, fill = FALSE,
        sv_box("R²", sprintf("%.0f%%", 100 * s$r.squared), p("share of variation explained"), tone = "blue"),
        sv_box("Adjusted R²", sprintf("%.0f%%", 100 * s$adj.r.squared),
               p("R² with a penalty per predictor: it only rises if a predictor earns its place"), tone = "gray")
      )
    })

    output$real_resid <- renderPlot({
      f <- real_fit()
      ggplot(data.frame(fv = fitted(f), r = resid(f)), aes(fv, r)) +
        geom_hline(yintercept = 0, colour = sv_col$ink2) +
        geom_point(colour = sv_col$primary, alpha = 0.6, size = 1.8) +
        geom_smooth(method = "loess", formula = y ~ x, se = FALSE, colour = sv_col$accent, linewidth = 0.9) +
        labs(x = "Fitted value", y = "Residual") + theme_statviz(11)
    })

    output$coef_sentences <- renderUI({
      df <- coef_df(); yv <- input$yvar
      others <- function(t) setdiff(df$term, t)
      unit <- if (isTRUE(input$standardize)) "one standard deviation" else "one unit"
      tags$ul(lapply(seq_len(nrow(df)), function(i) {
        r <- df[i, ]
        tags$li(
          strong(r$term), ": each ", unit, " higher ", r$term, " goes with ",
          strong(sprintf("%s %s", fmt(abs(r$est), 3), if (r$est >= 0) "higher" else "lower")), " ", yv,
          if (length(others(r$term))) paste0(", holding ", paste(others(r$term), collapse = ", "), " fixed"), ".",
          span(class = "text-muted",
               sprintf(" 95%% CI %s to %s, p %s.", fmt(r$lo, 3), fmt(r$hi, 3),
                       if (r$p < 0.001) "< 0.001" else paste("=", fmt(r$p, 3))),
               if (r$lo < 0 && r$hi > 0) "The interval includes 0: no clear evidence of an effect." else "")
        )
      }))
    })

    quiz_server(input, output, session, list(
      list(q = "In y = b0 + b1·x1 + b2·x2, the coefficient b1 tells you...",
           choices = c("the correlation between x1 and y",
                       "the change in y per unit of x1, holding x2 fixed",
                       "the change in y per unit of x1, ignoring x2",
                       "how much of y's variation x1 explains"),
           answer = 2,
           why = "Each coefficient is a 'holding the others constant' effect, as the added-variable plot shows."),
      list(q = "x1 has no real effect on y, but x1 is strongly correlated with x2, which does. A simple regression of y on x1 will probably show...",
           choices = c("a slope near zero", "a misleadingly strong slope", "a negative R²", "an error"),
           answer = 2,
           why = "x1 soaks up x2's effect: omitted-variable bias (confounding)."),
      list(q = "Strong correlation between predictors mainly causes...",
           choices = c("bad predictions", "biased coefficients", "large standard errors for the coefficients",
                       "R² above 1"),
           answer = 3,
           why = "Multicollinearity makes it hard to split credit between predictors, so the coefficients become unstable. Predictions can still be fine."),
      list(q = "Why look at adjusted R² instead of R²?",
           choices = c("R² never goes down when you add a predictor, even a useless one",
                       "Adjusted R² is always larger", "Adjusted R² ignores the intercept",
                       "They're the same thing"),
           answer = 1,
           why = "Adjusted R² charges a penalty for each predictor, so it only rises if the new predictor helps more than chance would.")
    ))
  })
}
