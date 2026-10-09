# Logistic regression module ---------------------------------------------------
#
# Tabs:
#   1. Why not a straight line? - linear vs logistic fit to 0/1 data
#   2. Odds & log-odds          - the scales, and how to read the coefficients
#   3. Fitting by likelihood    - fit by hand; log-likelihood landscape
#   4. From probabilities to decisions - threshold, confusion matrix, ROC
#   5. Real data
#   6. Quiz

logreg_ui <- function(id) {
  ns <- NS(id)

  layout_sidebar(
    sidebar = sidebar(
      width = 300,
      show_on(ns, c("why", "odds", "fit", "classify"),
        h6("Simulated data"),
        p(class = "small text-muted",
          "Students' hours of study (x) and whether they passed an exam (y = 1) or failed (y = 0).",
          "The chance of passing follows a true S-curve that you control."),
        sliderInput(ns("n"), "Number of students", 20, 500, 100, step = 10),
        sliderInput(ns("mid"), "True 50/50 point (hours)", 1, 9, 5, step = 0.5),
        sliderInput(ns("steep"), "True steepness of the curve", 0.2, 3, 1, step = 0.1),
        actionButton(ns("resample"), "Draw a new sample", icon = icon("shuffle"),
                     class = "btn-outline-secondary btn-sm w-100")
      ),
      show_on(ns, "odds",
        hr(),
        sliderInput(ns("p_demo"), "Try a probability", 0.01, 0.99, 0.75, step = 0.01)
      ),
      show_on(ns, "fit",
        hr(),
        h6("Your curve"),
        sliderInput(ns("my_b0"), "Intercept b0", -12, 4, -2, step = 0.1),
        sliderInput(ns("my_b1"), "Slope b1", -1, 4, 0.5, step = 0.05),
        actionButton(ns("reveal"), "Show the best (maximum-likelihood) curve", icon = icon("bullseye"),
                     class = "btn-primary w-100")
      ),
      show_on(ns, "classify",
        hr(),
        sliderInput(ns("thresh"), "Decision threshold", 0.05, 0.95, 0.5, step = 0.01),
        helpText("Predict 'pass' when the predicted probability is at least this.")
      ),
      show_on(ns, "real",
        selectInput(ns("dataset"), "Dataset",
                    setNames(names(logit_datasets), vapply(logit_datasets, `[[`, "", "name"))),
        uiOutput(ns("blurb")),
        selectizeInput(ns("xvars"), "Predictors", NULL, multiple = TRUE,
                       options = list(plugins = list("remove_button"))),
        selectInput(ns("effect_var"), "Show the effect of", NULL)
      ),
      show_on(ns, "quiz", p("Test yourself on the ideas from the other tabs."))
    ),

    navset_underline(
      id = ns("tab"),

      # ==== 1. Why ============================================================
      nav_panel(
        "1 · Why not a straight line?", value = "why",
        explain(
          p("Sometimes the outcome is yes/no: passed or failed, sick or healthy, bought or didn't.",
            "We code it as 1 or 0 and want to predict the ", strong("probability"), " of a 1."),
          p("A straight line (orange) is a poor fit: it runs below 0 and above 1, which no probability can do,",
            "and real probabilities usually level off at the extremes (studying 20 hours instead of 15",
            "barely changes your chances once you're almost certain to pass)."),
          p(strong("Logistic regression"), " fits an S-shaped curve (blue) that always stays between 0 and 1.",
            "The next tab shows the trick that makes this work.")
        ),
        card(full_screen = TRUE, card_header("Pass (1) or fail (0) against hours studied"),
             plotOutput(ns("why_plot"), height = "430px"),
             axes_note("hours studied (x).",
                       "outcome: 1 = passed, 0 = failed (points are nudged up/down slightly so they don't overlap).",
                       "Read the curves as 'probability of passing'. The orange line predicts impossible values",
                       "outside the shaded 0-1 zone."))
      ),

      # ==== 2. Odds ===========================================================
      nav_panel(
        "2 · Odds & log-odds", value = "odds",
        explain(withMathJax(
          p("Probabilities are stuck between 0 and 1, so we can't model them with a straight line directly.",
            "The fix is to stretch them onto a scale that runs from minus infinity to plus infinity, in two steps:"),
          math_words(
            "$$\\text{odds} = \\frac{p}{1-p}$$",
            "the chance it happens divided by the chance it doesn't. A probability of 0.75 is odds of",
            "0.75 / 0.25 = 3, often said '3 to 1'. Odds run from 0 to infinity.",
            symbols = list("\\(p\\)" = "probability of a 1 (e.g. passing)")
          ),
          math_words(
            "$$\\text{log-odds} = \\log\\Big(\\frac{p}{1-p}\\Big)$$",
            "taking the (natural) log stretches the odds over the whole number line: p = 0.5 becomes 0,",
            "likely events become positive, unlikely ones negative. This is also called the ", strong("logit"), "."
          ),
          p(strong("Logistic regression is a straight line on the log-odds scale:")),
          math_words(
            "$$\\log\\Big(\\frac{p}{1-p}\\Big) = b_0 + b_1\\,x$$",
            "each extra unit of x adds b1 to the log-odds. Undo the log and that means each extra unit of x ",
            strong("multiplies the odds by"), " \\(e^{b_1}\\), a number called the ", strong("odds ratio"), ".",
            symbols = list("\\(e^{b_1}\\)" = "'e to the power b1' (exp(b1) in R)")
          )
        )),
        layout_columns(
          col_widths = c(5, 7),
          div(
            uiOutput(ns("odds_boxes")),
            card(card_header("Probability vs log-odds"),
                 plotOutput(ns("logit_curve"), height = "260px"),
                 axes_note("log-odds.", "probability.",
                           "This S-curve converts any log-odds back into a probability. The orange dot is the",
                           "probability you chose in the sidebar."))
          ),
          card(full_screen = TRUE, card_header("The fitted model on both scales"),
               plotOutput(ns("two_scales"), height = "360px"),
               axes_note("hours studied.", "left: log-odds of passing; right: probability of passing.",
                         "Same model, two views. On the left it's a straight line; on the right the same line",
                         "becomes an S-curve."),
               uiOutput(ns("odds_interp")))
        )
      ),

      # ==== 3. Fitting ========================================================
      nav_panel(
        "3 · Fitting by likelihood", value = "fit",
        explain(withMathJax(
          p("There are no vertical residuals to square here (each y is just 0 or 1), so least squares isn't",
            "the natural choice. Instead we ask: ",
            em("if this curve were true, how likely would it be to see exactly the outcomes we saw?")),
          p("For each student, the curve gives a probability of passing. If they passed, the curve 'gets'",
            "that probability; if they failed, it gets 1 minus it. A good curve gives high probability to",
            "what actually happened. The gray segments show what each student's probability falls short of 1."),
          math_words(
            "$$\\log L = \\sum_i \\Big[ y_i \\log p_i + (1-y_i)\\log(1-p_i) \\Big]$$",
            "for each point, take the log of the probability the curve gave to what really happened,",
            "and add them up. This ", strong("log-likelihood"), " is always negative; closer to 0 is better.",
            "The best curve has the highest log-likelihood: the ", strong("maximum-likelihood"), " estimate.",
            "(Logs turn 'multiply all the probabilities' into 'add them up', which is easier to work with.)",
            symbols = list("\\(p_i\\)" = "the curve's probability of a 1 for point i",
                           "\\(y_i\\)" = "1 or 0: what happened")
          ),
          p("Unlike least squares, there's no formula for the answer. The computer climbs the hill in the",
            "landscape step by step until it reaches the top.")
        )),
        layout_columns(
          col_widths = c(7, 5),
          card(full_screen = TRUE, card_header("The data and your curve"),
               plotOutput(ns("fit_plot"), height = "430px"),
               axes_note("hours studied.", "probability of passing; points at 1 passed, at 0 failed.",
                         "Orange: your curve. Gray segments: how far each point's probability falls short of 1.",
                         "Green dashed: the maximum-likelihood curve, once revealed.")),
          div(
            uiOutput(ns("fit_stats")),
            card(card_header("The likelihood landscape"),
                 plotOutput(ns("ll_map"), height = "300px"),
                 axes_note("intercept b0.", "slope b1.",
                           "Darker = higher log-likelihood (better). Orange: your curve. Green: the peak."))
          )
        )
      ),

      # ==== 4. Classify =======================================================
      nav_panel(
        "4 · From probabilities to decisions", value = "classify",
        explain(
          p("A probability isn't a decision. To predict pass or fail, pick a ", strong("threshold"),
            ": predict 'pass' if the probability is at least that. Every threshold makes two kinds of mistake:"),
          tags$ul(
            tags$li(strong("False positive:"), " predicting pass for someone who failed."),
            tags$li(strong("False negative:"), " predicting fail for someone who passed.")
          ),
          p(strong("Sensitivity"), " is the share of real passes we catch; ", strong("specificity"),
            " is the share of real fails we correctly call fails. Lowering the threshold raises",
            "sensitivity but lowers specificity. The ", strong("ROC curve"), " shows that trade-off at every threshold."),
          p("The area under it (", strong("AUC"), ") summarizes how well the model separates the two groups:",
            "0.5 = no better than a coin flip, 1 = perfect.")
        ),
        layout_columns(
          col_widths = c(5, 3, 4),
          card(card_header("Predicted probability by actual outcome"),
               plotOutput(ns("class_strip"), height = "320px"),
               axes_note("the model's predicted probability of passing.", "what really happened.",
                         "Dashed line: your threshold. Points right of it are predicted 'pass'.")),
          card(card_header("Confusion matrix"),
               plotOutput(ns("confusion"), height = "260px"),
               uiOutput(ns("class_metrics"))),
          card(card_header("ROC curve"),
               plotOutput(ns("roc"), height = "320px"),
               axes_note("false positive rate (1 − specificity).", "sensitivity (true positive rate).",
                         "Each point on the curve is one threshold. Orange dot: your threshold."))
        )
      ),

      # ==== 5. Real data ======================================================
      nav_panel(
        "5 · Real data", value = "real",
        explain(p(
          "Pick a dataset and predictors. Coefficients are reported as odds ratios: the number the odds get",
          "multiplied by for each 1-unit increase in a predictor, holding the others fixed.",
          "An odds ratio of 1 means no effect; above 1, the outcome becomes more likely; below 1, less likely."
        )),
        layout_columns(
          col_widths = c(6, 6),
          card(card_header("Odds ratios with 95% confidence intervals"),
               plotOutput(ns("or_plot"), height = "320px"),
               axes_note("odds ratio per 1-unit increase (log scale, so ×2 and ÷2 look equally far from 1).",
                         "the predictors.",
                         "If a bar crosses the dashed line at 1, the data are consistent with no effect.")),
          card(card_header("Predicted probability"),
               plotOutput(ns("effect_plot"), height = "320px"),
               axes_note("the chosen predictor.", "predicted probability of the outcome.",
                         "Other predictors are held at their averages. Band: 95% confidence interval."))
        ),
        card(card_header("What each coefficient means"), uiOutput(ns("or_sentences")))
      ),

      nav_panel("Quiz", value = "quiz", quiz_panel(ns))
    )
  )
}

logreg_server <- function(id) {
  moduleServer(id, function(input, output, session) {

    seed <- reactiveVal(31)
    observeEvent(input$resample, seed(sample.int(1e6, 1)))
    true_b <- reactive(c(-input$mid * input$steep, input$steep))
    dat <- reactive(sim_logistic(input$n, true_b()[1], true_b()[2], seed()))
    fit <- reactive(glm(y ~ x, family = binomial, data = dat()))
    mle <- reactive(unname(coef(fit())))
    jit <- reactive({ set.seed(1); runif(nrow(dat()), -0.03, 0.03) })

    xs <- seq(0, 10, length.out = 200)

    # ---- 1. why -------------------------------------------------------------
    output$why_plot <- renderPlot({
      d <- dat(); lin <- coef(lm(y ~ x, data = d))
      curve <- data.frame(x = xs, p = plogis(mle()[1] + mle()[2] * xs))
      ggplot() +
        annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = 1, fill = sv_col$neutral, alpha = 0.6) +
        geom_hline(yintercept = c(0, 1), colour = sv_col$ink2, linetype = "dotted") +
        geom_abline(intercept = lin[1], slope = lin[2], colour = sv_col$accent, linewidth = 1.2) +
        geom_line(data = curve, aes(x, p), colour = sv_col$primary, linewidth = 1.4) +
        geom_point(data = d, aes(x, y + jit()), colour = sv_col$ink2, alpha = 0.5, size = 2) +
        annotate("label", x = 9.8, y = lin[1] + lin[2] * 9.8, label = "straight line", hjust = 1, vjust = -0.3,
                 colour = sv_col$accent, fill = sv_col$surface, label.size = 0, fontface = "bold") +
        annotate("label", x = 9.8, y = plogis(mle()[1] + mle()[2] * 9.8), label = "logistic curve",
                 hjust = 1, vjust = 1.3, colour = sv_col$primary, fill = sv_col$surface, label.size = 0,
                 fontface = "bold") +
        coord_cartesian(xlim = c(0, 10), ylim = c(-0.35, 1.35)) +
        scale_y_continuous(breaks = seq(0, 1, 0.25)) +
        labs(x = "Hours studied", y = "Pass (1) / fail (0)") + theme_statviz()
    })

    # ---- 2. odds ------------------------------------------------------------
    output$odds_boxes <- renderUI({
      p <- input$p_demo; o <- p / (1 - p)
      layout_column_wrap(
        width = 1 / 3, fill = FALSE,
        sv_box("Probability", fmt(p), tone = "blue"),
        sv_box("Odds", fmt(o), p(if (o >= 1) sprintf("%s to 1", fmt(o, 1)) else sprintf("1 to %s", fmt(1 / o, 1))),
               tone = "blue"),
        sv_box("Log-odds", fmt(log(o)), tone = "blue")
      )
    })

    output$logit_curve <- renderPlot({
      lo <- seq(-6, 6, length.out = 200)
      p <- input$p_demo
      ggplot(data.frame(lo, p = plogis(lo)), aes(lo, p)) +
        geom_hline(yintercept = c(0, 1), colour = sv_col$grid) +
        geom_line(colour = sv_col$primary, linewidth = 1.2) +
        annotate("segment", x = qlogis(p), xend = qlogis(p), y = 0, yend = p, colour = sv_col$accent, linetype = "22") +
        annotate("segment", x = -6, xend = qlogis(p), y = p, yend = p, colour = sv_col$accent, linetype = "22") +
        annotate("point", x = qlogis(p), y = p, colour = sv_col$accent, size = 4) +
        labs(x = "Log-odds", y = "Probability") + theme_statviz()
    })

    output$two_scales <- renderPlot({
      b <- mle()
      df <- rbind(data.frame(x = xs, v = b[1] + b[2] * xs, panel = "Log-odds scale: a straight line"),
                  data.frame(x = xs, v = plogis(b[1] + b[2] * xs), panel = "Probability scale: an S-curve"))
      df$panel <- factor(df$panel, levels = unique(df$panel))
      ggplot(df, aes(x, v)) +
        geom_line(colour = sv_col$primary, linewidth = 1.2) +
        geom_path(data = data.frame(x = c(4, 5, 5), v = b[1] + b[2] * c(4, 4, 5),
                                    panel = factor(levels(df$panel)[1], levels = levels(df$panel))),
                  colour = sv_col$accent, linewidth = 1) +
        geom_text(data = data.frame(x = 5.15, v = b[1] + b[2] * 4.5, lab = sprintf("+%s per hour", fmt(b[2])),
                                    panel = factor(levels(df$panel)[1], levels = levels(df$panel))),
                  aes(label = lab), hjust = 0, colour = sv_col$accent, size = 3.6) +
        facet_wrap(~panel, scales = "free_y") +
        labs(x = "Hours studied", y = NULL) +
        theme_statviz() + theme(strip.text = element_text(face = "bold", colour = sv_col$ink))
    })

    output$odds_interp <- renderUI({
      b <- mle(); or <- exp(b[2])
      div(class = "sv-math small mx-3",
          p(class = "mt-2 mb-1", strong("Reading this sample's fit:")),
          tags$ul(class = "mb-1",
            tags$li(sprintf("Slope b1 = %s: each extra hour adds %s to the log-odds of passing.", fmt(b[2]), fmt(b[2]))),
            tags$li(sprintf("Odds ratio e^b1 = %s: each extra hour multiplies the odds of passing by %s",
                            fmt(or), fmt(or)),
                    sprintf("(%s%.0f%%).", if (or >= 1) "+" else "", 100 * (or - 1))),
            tags$li(sprintf("50/50 point: −b0 / b1 = %s hours (true value %s).", fmt(-b[1] / b[2], 1), fmt(input$mid, 1)))
          ),
          p(class = "text-muted mb-1",
            "Note: the odds ratio is constant, but the change in ", em("probability"),
            " per hour is not. It's biggest near the 50/50 point and tiny at the flat ends of the S."))
    })

    # ---- 3. fitting ---------------------------------------------------------
    revealed <- reactiveVal(FALSE)
    observeEvent(list(input$n, input$mid, input$steep, input$resample), revealed(FALSE), ignoreInit = TRUE)
    observeEvent(input$reveal, {
      revealed(TRUE)
      updateSliderInput(session, "my_b0", value = round(max(-12, min(4, mle()[1])), 1))
      updateSliderInput(session, "my_b1", value = round(max(-1, min(4, mle()[2])) / 0.05) * 0.05)
    })
    loglik <- function(d, b0, b1) {
      eta <- b0 + b1 * d$x
      sum(d$y * eta - log1p(exp(eta)))
    }

    output$fit_plot <- renderPlot({
      d <- dat(); b0 <- input$my_b0; b1 <- input$my_b1
      d$p <- plogis(b0 + b1 * d$x)
      d$yj <- d$y + jit()
      g <- ggplot(d) +
        geom_segment(aes(x, yj, xend = x, yend = p), colour = sv_col$muted, linewidth = 0.4, alpha = 0.7)
      if (revealed()) {
        g <- g + geom_line(data = data.frame(x = xs, p = plogis(mle()[1] + mle()[2] * xs)), aes(x, p),
                           colour = sv_col$good, linewidth = 1, linetype = "22")
      }
      g + geom_line(data = data.frame(x = xs, p = plogis(b0 + b1 * xs)), aes(x, p),
                    colour = sv_col$accent, linewidth = 1.3) +
        geom_point(aes(x, yj), colour = sv_col$primary, alpha = 0.7, size = 2.2) +
        coord_cartesian(ylim = c(-0.05, 1.05)) +
        labs(x = "Hours studied", y = "Probability of passing") + theme_statviz()
    })

    output$fit_stats <- renderUI({
      d <- dat(); n <- nrow(d)
      mine <- loglik(d, input$my_b0, input$my_b1)
      best <- loglik(d, mle()[1], mle()[2])
      layout_column_wrap(
        width = 1 / 2, fill = FALSE,
        sv_box("Your log-likelihood", fmt(mine, 1),
               p(sprintf("typical probability given to what happened: %.0f%%", 100 * exp(mine / n))), tone = "orange"),
        sv_box(if (revealed()) "Maximum log-likelihood" else "Best possible",
               if (revealed()) fmt(best, 1) else "?",
               p(if (revealed()) sprintf("b0 = %s, b1 = %s (typical: %.0f%%)", fmt(mle()[1]), fmt(mle()[2]),
                                         100 * exp(best / n))
                 else "press the button to reveal"),
               tone = "green")
      )
    })

    output$ll_map <- renderPlot({
      d <- dat()
      g0 <- seq(-12, 4, length.out = 80); g1 <- seq(-1, 4, length.out = 80)
      grid <- expand.grid(b0 = g0, b1 = g1)
      eta <- outer(d$x, grid$b1) + matrix(grid$b0, nrow(d), nrow(grid), byrow = TRUE)
      grid$ll <- colSums(d$y * eta - log1p(exp(eta)))
      best <- loglik(d, mle()[1], mle()[2])
      grid$v <- pmax(grid$ll, best - 60)         # keep the colour scale focused near the peak
      g <- ggplot(grid, aes(b0, b1)) +
        geom_raster(aes(fill = v), interpolate = TRUE) +
        geom_contour(aes(z = v), colour = "#ffffff", alpha = 0.5, bins = 12, linewidth = 0.3) +
        scale_fill_gradientn(colours = seq_blues, name = "log L") +
        annotate("point", x = input$my_b0, y = input$my_b1, colour = sv_col$surface, size = 5.5) +
        annotate("point", x = input$my_b0, y = input$my_b1, colour = sv_col$accent, size = 4)
      if (revealed()) {
        g <- g + annotate("point", x = mle()[1], y = mle()[2], colour = sv_col$surface, size = 5.5) +
          annotate("point", x = mle()[1], y = mle()[2], colour = sv_col$good, size = 4, shape = 18)
      }
      g + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
        labs(x = "Intercept b0", y = "Slope b1") + theme_statviz() + theme(legend.position = "right")
    })

    # ---- 4. classify --------------------------------------------------------
    probs <- reactive(fitted(fit()))

    output$class_strip <- renderPlot({
      d <- dat(); d$p <- probs()
      d$actual <- factor(ifelse(d$y == 1, "Passed", "Failed"), levels = c("Failed", "Passed"))
      set.seed(2)
      ggplot(d, aes(p, actual, colour = actual)) +
        geom_vline(xintercept = input$thresh, colour = sv_col$ink, linetype = "22", linewidth = 0.9) +
        geom_jitter(height = 0.2, width = 0, alpha = 0.7, size = 2.2) +
        scale_colour_manual(values = c(Failed = sv_col$accent, Passed = sv_col$primary), guide = "none") +
        scale_x_continuous(limits = c(0, 1)) +
        labs(x = "Predicted probability of passing", y = NULL) + theme_statviz()
    })

    cm <- reactive({
      pred <- probs() >= input$thresh; y <- dat()$y == 1
      c(TP = sum(pred & y), FP = sum(pred & !y), FN = sum(!pred & y), TN = sum(!pred & !y))
    })

    output$confusion <- renderPlot({
      m <- cm()
      df <- data.frame(pred = c("Pass", "Pass", "Fail", "Fail"), actual = c("Passed", "Failed", "Passed", "Failed"),
                       n = c(m["TP"], m["FP"], m["FN"], m["TN"]),
                       lab = c("true positive", "false positive", "false negative", "true negative"),
                       ok = c(TRUE, FALSE, FALSE, TRUE))
      ggplot(df, aes(pred, actual)) +
        geom_tile(aes(fill = ok), colour = sv_col$surface, linewidth = 2) +
        geom_text(aes(label = n), size = 7, fontface = "bold", colour = sv_col$ink, vjust = -0.1) +
        geom_text(aes(label = lab), size = 3.2, colour = sv_col$ink2, vjust = 2.2) +
        scale_fill_manual(values = c(`TRUE` = "#e6f6ef", `FALSE` = "#fdeee7"), guide = "none") +
        labs(x = "Predicted", y = "Actual") + theme_statviz() + theme(panel.grid = element_blank())
    })

    output$class_metrics <- renderUI({
      m <- cm(); n <- sum(m)
      tags$table(
        class = "table table-sm small mb-0",
        tags$tr(tags$td("Accuracy"), tags$td(strong(sprintf("%.0f%%", 100 * (m["TP"] + m["TN"]) / n)))),
        tags$tr(tags$td("Sensitivity"), tags$td(sprintf("%.0f%%", 100 * m["TP"] / max(1, m["TP"] + m["FN"])))),
        tags$tr(tags$td("Specificity"), tags$td(sprintf("%.0f%%", 100 * m["TN"] / max(1, m["TN"] + m["FP"]))))
      )
    })

    output$roc <- renderPlot({
      p <- probs(); y <- dat()$y
      th <- c(Inf, sort(unique(p), decreasing = TRUE), -Inf)
      roc <- data.frame(fpr = vapply(th, function(t) mean(p[y == 0] >= t), 0),
                        tpr = vapply(th, function(t) mean(p[y == 1] >= t), 0))
      auc <- sum(diff(roc$fpr) * (head(roc$tpr, -1) + tail(roc$tpr, -1)) / 2)
      cur <- c(mean(p[y == 0] >= input$thresh), mean(p[y == 1] >= input$thresh))
      ggplot(roc, aes(fpr, tpr)) +
        geom_abline(colour = sv_col$muted, linetype = "22") +
        geom_area(fill = "#cde2fb", alpha = 0.6) +
        geom_step(colour = sv_col$primary, linewidth = 1.1, direction = "vh") +
        annotate("point", x = cur[1], y = cur[2], colour = sv_col$surface, size = 5.5) +
        annotate("point", x = cur[1], y = cur[2], colour = sv_col$accent, size = 4) +
        annotate("label", x = 0.6, y = 0.25, label = sprintf("AUC = %s", fmt(auc)),
                 colour = sv_col$ink, fill = sv_col$surface, label.size = 0, fontface = "bold") +
        coord_equal() +
        labs(x = "False positive rate", y = "Sensitivity") + theme_statviz()
    })

    # ---- 5. real data -------------------------------------------------------
    ds <- reactive(logit_datasets[[input$dataset]])
    output$blurb <- renderUI(p(class = "small text-muted", ds()$blurb,
                               br(), "Outcome: ", code(ds()$y)))
    observeEvent(ds(), {
      updateSelectizeInput(session, "xvars", choices = setdiff(names(ds()$data), ds()$y), selected = ds()$x)
    })
    observeEvent(input$xvars, {
      updateSelectInput(session, "effect_var", choices = input$xvars,
                        selected = if (isTRUE(input$effect_var %in% input$xvars)) input$effect_var else input$xvars[1])
    }, ignoreNULL = FALSE)

    real_fit <- reactive({
      d <- ds()$data
      xv <- intersect(input$xvars, names(d))
      validate(need(length(xv) >= 1, "Pick at least one predictor."))
      glm(reformulate(xv, response = ds()$y), family = binomial, data = d)
    })

    or_df <- reactive({
      f <- real_fit()
      s <- summary(f)$coefficients[-1, , drop = FALSE]
      ci <- confint.default(f)[-1, , drop = FALSE]
      data.frame(term = rownames(s), or = exp(s[, 1]), lo = exp(ci[, 1]), hi = exp(ci[, 2]), p = s[, 4])
    })

    output$or_plot <- renderPlot({
      df <- or_df(); df$term <- factor(df$term, levels = rev(df$term))
      ggplot(df, aes(or, term)) +
        geom_vline(xintercept = 1, colour = sv_col$ink2, linetype = "22") +
        geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.2, colour = sv_col$primary, linewidth = 0.9) +
        geom_point(colour = sv_col$primary, size = 3.5) +
        scale_x_log10() +
        labs(x = "Odds ratio (log scale)", y = NULL) + theme_statviz()
    })

    output$effect_plot <- renderPlot({
      f <- real_fit(); d <- ds()$data
      v <- input$effect_var
      req(v %in% input$xvars)
      xv <- intersect(input$xvars, names(d))
      nd <- as.data.frame(lapply(d[xv], function(col) rep(mean(col, na.rm = TRUE), 100)))
      nd[[v]] <- seq(min(d[[v]], na.rm = TRUE), max(d[[v]], na.rm = TRUE), length.out = 100)
      pr <- predict(f, nd, type = "link", se.fit = TRUE)
      nd$p <- plogis(pr$fit); nd$lo <- plogis(pr$fit - 1.96 * pr$se.fit); nd$hi <- plogis(pr$fit + 1.96 * pr$se.fit)
      ggplot(nd, aes(.data[[v]], p)) +
        geom_ribbon(aes(ymin = lo, ymax = hi), fill = sv_col$primary, alpha = 0.2) +
        geom_line(colour = sv_col$primary, linewidth = 1.2) +
        geom_rug(data = d, aes(x = .data[[v]], y = NULL), sides = "b", alpha = 0.3, colour = sv_col$ink2) +
        scale_y_continuous(limits = c(0, 1)) +
        labs(x = v, y = sprintf("Probability that %s = 1", ds()$y)) + theme_statviz()
    })

    output$or_sentences <- renderUI({
      df <- or_df()
      tags$ul(lapply(seq_len(nrow(df)), function(i) {
        r <- df[i, ]
        others <- setdiff(df$term, r$term)
        tags$li(
          strong(r$term), sprintf(": each 1-unit increase multiplies the odds of %s = 1 by ", ds()$y),
          strong(fmt(r$or, 3)),
          sprintf(" (%s%.0f%%)", if (r$or >= 1) "+" else "", 100 * (r$or - 1)),
          if (length(others)) paste0(", holding ", paste(others, collapse = ", "), " fixed"), ".",
          span(class = "text-muted",
               sprintf(" 95%% CI %s to %s, p %s.", fmt(r$lo, 3), fmt(r$hi, 3),
                       if (r$p < 0.001) "< 0.001" else paste("=", fmt(r$p, 3))),
               if (r$lo < 1 && r$hi > 1) "The interval includes 1: no clear evidence of an effect." else "")
        )
      }))
    })

    quiz_server(input, output, session, list(
      list(q = "Why not fit an ordinary straight line to a 0/1 outcome?",
           choices = c("It's too slow to compute", "It can predict probabilities below 0 or above 1",
                       "It can't use more than one predictor", "It always has R² = 0"),
           answer = 2,
           why = "A straight line keeps going forever; probabilities must stay between 0 and 1."),
      list(q = "A probability of 0.8 corresponds to odds of...",
           choices = c("0.8", "1.25", "4", "8"),
           answer = 3,
           why = "Odds = p / (1 − p) = 0.8 / 0.2 = 4, or '4 to 1'."),
      list(q = "In logistic regression, an odds ratio of 1.5 for x means...",
           choices = c("each unit of x adds 1.5 to the probability",
                       "each unit of x multiplies the odds by 1.5",
                       "the probability is 1.5 times higher overall",
                       "x explains 150% of the outcome"),
           answer = 2,
           why = "exp(b1) is the factor the odds get multiplied by per 1-unit increase in x."),
      list(q = "Logistic regression coefficients are found by...",
           choices = c("minimizing the sum of squared residuals",
                       "maximizing the likelihood of the observed outcomes",
                       "averaging the outcomes",
                       "a closed-form formula like the least-squares slope"),
           answer = 2,
           why = "Maximum likelihood: find the curve that makes what actually happened most probable. There's no closed form, so the computer climbs the likelihood hill step by step."),
      list(q = "Lowering the decision threshold from 0.5 to 0.3 will usually...",
           choices = c("raise sensitivity and lower specificity",
                       "lower sensitivity and raise specificity",
                       "raise both", "change the AUC"),
           answer = 1,
           why = "More cases get called positive: you catch more real positives but also make more false alarms. The AUC summarizes all thresholds, so it doesn't change.")
    ))
  })
}
