# Mixed-effects linear regression module ---------------------------------------
#
# Tabs:
#   1. Grouped data       - complete, no, and partial pooling
#   2. Shrinkage          - how the mixed model pulls group estimates together
#   3. Variance & honesty - variance components, ICC, standard errors
#   4. Real data: sleep study
#   5. Quiz

mixed_ui <- function(id) {
  ns <- NS(id)

  layout_sidebar(
    sidebar = sidebar(
      width = 300,
      show_on(ns, c("pool", "shrink", "variance"),
        h6("Simulated schools"),
        p(class = "small text-muted",
          "Students' test scores (y) against hours of homework (x), for students grouped in schools.",
          "Each school has its own baseline and its own homework effect."),
        sliderInput(ns("groups"), "Number of schools", 4, 30, 12, step = 1),
        sliderInput(ns("per_group"), "Students per school (average)", 2, 30, 8, step = 1),
        checkboxInput(ns("unbalanced"), "Schools differ in size", TRUE),
        sliderInput(ns("sd_int"), "How much school baselines differ (SD)", 0, 15, 8, step = 0.5),
        sliderInput(ns("sd_slope"), "How much the homework effect differs (SD)", 0, 3, 1, step = 0.1),
        sliderInput(ns("sd_res"), "Student-level noise (SD)", 1, 15, 6, step = 0.5),
        checkboxInput(ns("rand_slope"), "Let slopes vary by school in the model", TRUE),
        actionButton(ns("resample"), "Draw a new sample", icon = icon("shuffle"),
                     class = "btn-outline-secondary btn-sm w-100")
      ),
      show_on(ns, "sleep",
        radioButtons(ns("sleep_model"), "Model", c(
          "Each person has their own baseline" = "int",
          "Own baseline and own daily change" = "slope"
        ), selected = "slope")
      ),
      show_on(ns, "quiz", p("Test yourself on the ideas from the other tabs."))
    ),

    navset_underline(
      id = ns("tab"),

      # ==== 1. Pooling =========================================================
      nav_panel(
        "1 · Grouped data", value = "pool",
        explain(withMathJax(
          p("Data often come in groups: students in schools, patients in hospitals, repeated measurements",
            "on the same person. Observations in the same group tend to be more alike than observations",
            "from different groups, which breaks the 'independence' assumption of ordinary regression.",
            "There are three ways to handle it:"),
          tags$ul(
            tags$li(strong("Complete pooling"), " (gray dashed): ignore the groups and fit one line to everyone.",
                    "Simple, but it misses real differences between schools and is overconfident (tab 3)."),
            tags$li(strong("No pooling"), " (orange): fit a separate line for each school.",
                    "Respects differences, but a school with only 2 or 3 students gets a wild, unreliable line."),
            tags$li(strong("Partial pooling"), " (blue): a ", strong("mixed-effects model"),
                    ". Each school gets its own line, but the lines are pulled toward the overall average,",
                    "more strongly when a school has little data. It's a compromise that uses information from all",
                    "schools to improve each school's estimate.")
          ),
          math_words(
            "$$y_{ij} = (\\beta_0 + u_{0j}) + (\\beta_1 + u_{1j})\\,x_{ij} + \\varepsilon_{ij}$$",
            "student i in school j: the overall intercept and slope (the ", strong("fixed effects"), ", β), plus",
            "a nudge for school j (the ", strong("random effects"), ", u), plus that student's own noise (ε).",
            symbols = list("\\(\\beta_0, \\beta_1\\)" = "average intercept and slope across all schools",
                           "\\(u_{0j}, u_{1j}\\)" = "how school j differs from average",
                           "\\(\\varepsilon_{ij}\\)" = "student-level noise")
          ),
          math_words(
            "$$u_{0j} \\sim N(0,\\, \\tau_0^2) \\qquad u_{1j} \\sim N(0,\\, \\tau_1^2) \\qquad \\varepsilon_{ij} \\sim N(0,\\, \\sigma^2)$$",
            "the school nudges are assumed to come from a bell curve centered on 0. Schools differ, but",
            "most are near average and extreme schools are rare. That assumption is what lets a small",
            "school 'borrow strength' from the others.",
            symbols = list("\\(\\sim N(0, \\tau^2)\\)" = "'drawn from a bell curve with average 0 and spread τ'",
                           "\\(\\tau_0, \\tau_1\\)" = "how much schools vary (between-group SDs)",
                           "\\(\\sigma\\)" = "how much students vary within a school")
          )
        )),
        card(full_screen = TRUE,
             card_header("Each school's data and three fitted lines"),
             plotOutput(ns("pool_plot"), height = "620px"),
             axes_note("hours of homework.", "test score.",
                       "One panel per school (n = number of students). Gray dashed: complete pooling (the same in every",
                       "panel). Orange: no pooling. Blue: partial pooling (mixed model). Look at the smallest schools:",
                       "the orange line can be extreme, while the blue line stays sensible."))
      ),

      # ==== 2. Shrinkage =======================================================
      nav_panel(
        "2 · Shrinkage", value = "shrink",
        explain(
          p("The mixed model's school estimates sit between each school's own estimate and the overall average.",
            "This pull toward the middle is called ", strong("shrinkage"), "."),
          p("How far each school is pulled depends on how much we trust its own data. A big school with",
            "consistent students barely moves. A school with 2 students is pulled most of the way to the",
            "average, because its own estimate is mostly noise. On average, shrunken estimates are closer to",
            "the truth than the no-pooling ones. That's the payoff for 'borrowing strength'."),
          p("Try: make the schools very different (raise 'How much school baselines differ') and shrinkage",
            "weakens, because the average is less informative about any one school. Make them similar and",
            "everything collapses toward the average.")
        ),
        layout_columns(
          col_widths = c(7, 5),
          card(full_screen = TRUE, card_header("Each school's estimates: own data vs mixed model"),
               plotOutput(ns("shrink_plot"), height = "460px"),
               axes_note("predicted score at 5 hours of homework (the school's typical level).",
                         "homework effect: points gained per extra hour (the school's slope).",
                         "Orange hollow: no-pooling estimate. Blue: mixed-model estimate. Arrows show the pull.",
                         "Green cross: the overall average. Bigger circles = more students.")),
          div(
            card(card_header("How far each school moved"),
                 plotOutput(ns("shrink_size"), height = "260px"),
                 axes_note("number of students in the school.",
                           "distance between the school's own estimate and its mixed-model estimate (predicted score at 5 hours).",
                           "Small schools get pulled the most.")),
            uiOutput(ns("shrink_accuracy"))
          )
        )
      ),

      # ==== 3. Variance ========================================================
      nav_panel(
        "3 · Variance & honesty", value = "variance",
        explain(withMathJax(
          p("A mixed model splits the variation in scores into two parts: differences ", em("between"),
            " schools and differences between students ", em("within"), " a school."),
          math_words(
            "$$\\text{ICC} = \\frac{\\tau_0^2}{\\tau_0^2 + \\sigma^2}$$",
            "the ", strong("intraclass correlation"), " is the share of the total variation that comes from",
            "differences between schools. It's also how similar two students from the same school are",
            "expected to be. ICC = 0 means school doesn't matter; ICC = 1 means everyone in a school is identical.",
            symbols = list("\\(\\tau_0^2\\)" = "between-school variance (baselines)",
                           "\\(\\sigma^2\\)" = "within-school variance")
          ),
          p(strong("Why it matters for honesty:"), "complete pooling treats every student as an independent",
            "piece of evidence. But students in the same school share their school's quirks, so 100 students",
            "in 10 schools carry less information than 100 unrelated students. Ignoring that makes standard",
            "errors too small and confidence intervals too narrow, so you become overconfident. Use the",
            "simulation below to see how often each method's 95% interval really contains the truth.")
        )),
        layout_columns(
          col_widths = c(6, 6),
          card(card_header("Where the variation comes from"),
               plotOutput(ns("var_plot"), height = "220px"),
               axes_note("share of the total variance (random-intercept model).", "",
                         "Green: between schools. Gray: between students within a school."),
               uiOutput(ns("var_table"))),
          card(card_header("Estimating the homework effect"),
               uiOutput(ns("se_table")),
               actionButton(ns("run_cover"), "Run 100 simulated studies", icon = icon("play"),
                            class = "btn-outline-primary btn-sm", width = "260px"),
               uiOutput(ns("cover_result")))
        )
      ),

      # ==== 4. Sleep study =====================================================
      nav_panel(
        "4 · Real data: sleep study", value = "sleep",
        explain(p(
          "In this classic study, 18 people were limited to 3 hours of sleep per night, and their reaction",
          "times were measured every day for 10 days (days 0-1 were adaptation and baseline).",
          "Each person contributes 10 measurements, so the measurements are grouped by person: a perfect",
          "fit for a mixed model. (Data: ", code("lme4::sleepstudy"), ", Belenky et al. 2003.)"
        )),
        layout_columns(
          col_widths = c(8, 4),
          card(full_screen = TRUE, card_header("Each person's reaction times"),
               plotOutput(ns("sleep_plot"), height = "560px"),
               axes_note("days of sleep restriction.", "average reaction time (milliseconds; higher = slower).",
                         "One panel per person. Blue: that person's mixed-model line. Gray dashed: the average",
                         "person (fixed effects only).")),
          div(
            card(card_header("What the model says"), uiOutput(ns("sleep_summary"))),
            card(card_header("How people differ from average"),
                 plotOutput(ns("sleep_ranef"), height = "380px"),
                 axes_note("each person's deviation from the average (with a rough 95% interval).",
                           "people, sorted.",
                           "Left: baseline reaction time (ms). Right: extra ms per day. Intervals crossing 0",
                           "mean that person isn't clearly different from average."))
          )
        )
      ),

      nav_panel("Quiz", value = "quiz", quiz_panel(ns))
    )
  )
}

mixed_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    B0 <- 50; B1 <- 2

    quiet_lmer <- function(formula, data) {
      suppressMessages(suppressWarnings(lme4::lmer(formula, data = data, REML = TRUE)))
    }

    seed <- reactiveVal(41)
    observeEvent(input$resample, seed(sample.int(1e6, 1)))
    dat <- reactive(sim_mixed(input$groups, input$per_group, input$unbalanced, B0, B1,
                              input$sd_int, input$sd_slope, input$sd_res, seed()))

    model_formula <- reactive(if (isTRUE(input$rand_slope)) y ~ x + (x | group) else y ~ x + (1 | group))
    fit_mixed <- reactive(quiet_lmer(model_formula(), dat()))
    fit_pool <- reactive(lm(y ~ x, data = dat()))

    # per-group intercept/slope under each approach
    estimates <- reactive({
      d <- dat()
      grp <- levels(d$group)
      nop <- t(vapply(grp, function(g) {
        dg <- d[d$group == g, ]
        if (nrow(dg) >= 2 && sd(dg$x) > 0) coef(lm(y ~ x, data = dg)) else c(mean(dg$y) - B1 * mean(dg$x), NA)
      }, numeric(2)))
      pp <- as.matrix(coef(fit_mixed())$group[grp, c("(Intercept)", "x")])
      data.frame(group = grp, n = as.vector(table(d$group)[grp]),
                 a_no = nop[, 1], b_no = nop[, 2], a_pp = pp[, 1], b_pp = pp[, 2])
    })

    # ---- 1. pooling ---------------------------------------------------------
    output$pool_plot <- renderPlot({
      d <- dat(); e <- estimates(); cp <- coef(fit_pool())
      lab <- setNames(sprintf("%s (n = %d)", e$group, e$n), e$group)
      d$panel <- factor(lab[as.character(d$group)], levels = lab)
      e$panel <- factor(lab[e$group], levels = lab)
      ggplot(d, aes(x, y)) +
        geom_abline(intercept = cp[1], slope = cp[2], colour = sv_col$muted, linetype = "22", linewidth = 0.7) +
        geom_abline(data = e[!is.na(e$b_no), ], aes(intercept = a_no, slope = b_no),
                    colour = sv_col$accent, linewidth = 0.9) +
        geom_abline(data = e, aes(intercept = a_pp, slope = b_pp), colour = sv_col$primary, linewidth = 1.1) +
        geom_point(colour = sv_col$ink2, alpha = 0.6, size = 1.6) +
        facet_wrap(~panel, ncol = 6) +
        coord_cartesian(xlim = c(0, 10), ylim = range(d$y) + c(-5, 5)) +
        labs(x = "Hours of homework", y = "Test score") +
        theme_statviz(11) + theme(strip.text = element_text(colour = sv_col$ink, size = 9))
    })

    # ---- 2. shrinkage -------------------------------------------------------
    at5 <- reactive({
      e <- estimates()
      transform(e, lvl_no = a_no + 5 * ifelse(is.na(b_no), B1, b_no), lvl_pp = a_pp + 5 * b_pp,
                slope_no = ifelse(is.na(b_no), b_pp, b_no))
    })

    output$shrink_plot <- renderPlot({
      e <- at5(); fe <- lme4::fixef(fit_mixed())
      ggplot(e) +
        geom_segment(aes(lvl_no, slope_no, xend = lvl_pp, yend = b_pp), colour = sv_col$muted,
                     arrow = arrow(length = unit(0.18, "cm")), linewidth = 0.5) +
        geom_point(aes(lvl_no, slope_no, size = n), shape = 21, colour = sv_col$accent, fill = NA, stroke = 1.1) +
        geom_point(aes(lvl_pp, b_pp, size = n), colour = sv_col$primary, alpha = 0.85) +
        annotate("point", x = fe[1] + 5 * fe[2], y = fe[2], shape = 3, size = 7, stroke = 2, colour = sv_col$good) +
        scale_size_area(max_size = 7, name = "students") +
        labs(x = "Predicted score at 5 hours", y = "Points per extra hour") +
        theme_statviz() + theme(legend.position = "right")
    })

    output$shrink_size <- renderPlot({
      e <- at5()
      ggplot(e, aes(n, abs(lvl_no - lvl_pp))) +
        geom_point(colour = sv_col$primary, size = 2.6, alpha = 0.8) +
        labs(x = "Students in the school", y = "Amount of shrinkage") + theme_statviz(11)
    })

    output$shrink_accuracy <- renderUI({
      # compare each approach's school-level predictions with the true school lines
      tr <- attr(dat(), "truth")
      truth <- B0 + tr$u0 + 5 * (B1 + tr$u1)
      e <- at5()
      err_no <- sqrt(mean((e$lvl_no - truth)^2)); err_pp <- sqrt(mean((e$lvl_pp - truth)^2))
      sv_box("Typical error in a school's predicted score", fmt(err_pp, 1),
             p(sprintf("mixed model, vs %s using each school's data alone", fmt(err_no, 1)),
               br(), "(we know the truth because we simulated it)"),
             tone = if (err_pp <= err_no) "green" else "orange")
    })

    # ---- 3. variance --------------------------------------------------------
    fit_int <- reactive(quiet_lmer(y ~ x + (1 | group), dat()))

    output$var_plot <- renderPlot({
      vc <- as.data.frame(lme4::VarCorr(fit_int()))
      tau2 <- vc$vcov[vc$grp == "group"]; s2 <- vc$vcov[vc$grp == "Residual"]
      df <- data.frame(part = factor(c("Between schools", "Within schools"),
                                     levels = c("Within schools", "Between schools")),
                       v = c(tau2, s2) / (tau2 + s2))
      ggplot(df, aes(v, "", fill = part)) +
        geom_col(width = 0.5, colour = sv_col$surface, linewidth = 1) +
        geom_text(aes(label = sprintf("%s\n%.0f%%", part, 100 * v)), position = position_stack(vjust = 0.5),
                  colour = sv_col$ink, size = 3.8, fontface = "bold") +
        scale_fill_manual(values = c(`Between schools` = "#bfe9d6", `Within schools` = sv_col$neutral), guide = "none") +
        scale_x_continuous(labels = function(x) paste0(100 * x, "%")) +
        labs(x = NULL, y = NULL) + theme_statviz() + theme(panel.grid = element_blank())
    })

    output$var_table <- renderUI({
      vc <- as.data.frame(lme4::VarCorr(fit_mixed()))
      sdv <- function(var) { v <- vc$sdcor[vc$grp == var[1] & vc$var1 %in% var[2] & is.na(vc$var2)]; if (length(v)) v else NA }
      tau0 <- sdv(c("group", "(Intercept)")); tau1 <- sdv(c("group", "x"))
      sig <- vc$sdcor[vc$grp == "Residual"]
      vi <- as.data.frame(lme4::VarCorr(fit_int()))
      icc <- vi$vcov[vi$grp == "group"] / sum(vi$vcov)
      true_icc <- input$sd_int^2 / (input$sd_int^2 + input$sd_res^2)
      tagList(
        tags$table(
          class = "table table-sm small mt-2",
          tags$thead(tags$tr(tags$th(""), tags$th("Estimated"), tags$th("True"))),
          tags$tbody(
            tags$tr(tags$td("SD of school baselines (τ0)"), tags$td(fmt(tau0)), tags$td(fmt(input$sd_int))),
            tags$tr(tags$td("SD of school slopes (τ1)"),
                    tags$td(if (isTRUE(input$rand_slope)) fmt(tau1) else "not in model"), tags$td(fmt(input$sd_slope))),
            tags$tr(tags$td("SD of students within school (σ)"), tags$td(fmt(sig)), tags$td(fmt(input$sd_res))),
            tags$tr(tags$td(strong("ICC")), tags$td(strong(fmt(icc))), tags$td(fmt(true_icc)))
          )
        ),
        p(class = "small", sprintf("In words: about %.0f%% of the variation in scores is due to which school a student attends.",
                                   100 * icc))
      )
    })

    output$se_table <- renderUI({
      cp <- summary(fit_pool())$coefficients["x", ]
      mm <- summary(fit_mixed())$coefficients["x", ]
      row <- function(name, est, se) tags$tr(tags$td(name), tags$td(fmt(est)), tags$td(fmt(se, 3)),
                                           tags$td(sprintf("%s to %s", fmt(est - 1.96 * se), fmt(est + 1.96 * se))))
      tagList(
        tags$table(
          class = "table table-sm small",
          tags$thead(tags$tr(tags$th("Approach"), tags$th("Slope"), tags$th("SE"), tags$th("95% interval"))),
          tags$tbody(row("Complete pooling (lm)", cp[1], cp[2]),
                     row("Mixed model (lmer)", mm[1], mm[2]))
        ),
        p(class = "small text-muted", sprintf("True average slope: %s points per hour.", fmt(B1)),
          "The difference in SE is largest when slopes really vary between schools and the model allows it.")
      )
    })

    cover <- reactiveVal(NULL)
    observeEvent(list(input$groups, input$per_group, input$unbalanced, input$sd_int, input$sd_slope,
                      input$sd_res, input$rand_slope), cover(NULL), ignoreInit = TRUE)
    observeEvent(input$run_cover, {
      hits <- matrix(NA, 100, 2)
      withProgress(message = "Simulating studies", value = 0, {
        for (i in 1:100) {
          d <- sim_mixed(input$groups, input$per_group, input$unbalanced, B0, B1,
                         input$sd_int, input$sd_slope, input$sd_res, 1000 + i)
          a <- summary(lm(y ~ x, data = d))$coefficients["x", 1:2]
          b <- tryCatch(summary(quiet_lmer(model_formula(), d))$coefficients["x", 1:2], error = function(e) c(NA, NA))
          hits[i, ] <- c(abs(a[1] - B1) <= qt(0.975, nrow(d) - 2) * a[2], abs(b[1] - B1) <= 1.96 * b[2])
          incProgress(0.01)
        }
      })
      cover(colMeans(hits, na.rm = TRUE))
    })

    output$cover_result <- renderUI({
      cv <- cover()
      if (is.null(cv)) return(p(class = "small text-muted mt-2",
                                "Each simulated study draws new schools and students, fits both models,",
                                "and checks whether the 95% interval contains the true slope (2)."))
      layout_column_wrap(
        width = 1 / 2, fill = FALSE, class = "mt-2",
        sv_box("Complete pooling", sprintf("%.0f%%", 100 * cv[1]), p("of intervals contained the truth"),
               tone = if (cv[1] < 0.9) "orange" else "gray"),
        sv_box("Mixed model", sprintf("%.0f%%", 100 * cv[2]), p("should be close to 95%"),
               tone = if (cv[2] >= 0.9) "green" else "orange")
      )
    })

    # ---- 4. sleep study -----------------------------------------------------
    sleep_fit <- reactive({
      f <- if (input$sleep_model == "slope") Reaction ~ Days + (Days | Subject) else Reaction ~ Days + (1 | Subject)
      quiet_lmer(f, lme4::sleepstudy)
    })

    output$sleep_plot <- renderPlot({
      d <- lme4::sleepstudy; f <- sleep_fit()
      cf <- coef(f)$Subject; cf$Subject <- rownames(cf)
      fe <- lme4::fixef(f)
      ggplot(d, aes(Days, Reaction)) +
        geom_abline(intercept = fe[1], slope = fe[2], colour = sv_col$muted, linetype = "22", linewidth = 0.7) +
        geom_abline(data = cf, aes(intercept = `(Intercept)`, slope = Days), colour = sv_col$primary, linewidth = 1) +
        geom_point(colour = sv_col$ink2, size = 1.6, alpha = 0.8) +
        facet_wrap(~Subject, ncol = 6) +
        scale_x_continuous(breaks = c(0, 3, 6, 9)) +
        labs(x = "Days of sleep restriction", y = "Reaction time (ms)") +
        theme_statviz(11) + theme(strip.text = element_text(colour = sv_col$ink, size = 9))
    })

    output$sleep_summary <- renderUI({
      f <- sleep_fit()
      co <- summary(f)$coefficients
      vc <- as.data.frame(lme4::VarCorr(f))
      sd_int <- vc$sdcor[vc$grp == "Subject" & vc$var1 == "(Intercept)" & is.na(vc$var2)]
      sd_day <- vc$sdcor[vc$grp == "Subject" & vc$var1 == "Days" & is.na(vc$var2)]
      sig <- vc$sdcor[vc$grp == "Residual"]
      tagList(
        tags$ul(
          class = "small",
          tags$li(strong("Average baseline: "), sprintf("%s ms on day 0 (SE %s).", fmt(co[1, 1], 0), fmt(co[1, 2], 1))),
          tags$li(strong("Average effect of sleep loss: "),
                  sprintf("+%s ms per day (95%% CI roughly %s to %s).", fmt(co[2, 1], 1),
                          fmt(co[2, 1] - 1.96 * co[2, 2], 1), fmt(co[2, 1] + 1.96 * co[2, 2], 1))),
          tags$li(strong("People differ: "),
                  sprintf("baselines vary by about ±%s ms (SD)", fmt(sd_int, 0)),
                  if (length(sd_day)) sprintf(", and daily increases by about ±%s ms.", fmt(sd_day, 1)) else "."),
          tags$li(strong("Day-to-day noise: "), sprintf("±%s ms (SD) around each person's own line.", fmt(sig, 0)))
        ),
        if (input$sleep_model == "int")
          p(class = "small text-muted", "Every person's line has the same slope in this model. Compare the panels",
            "with the other model: several people clearly worsen faster or slower than average.")
      )
    })

    output$sleep_ranef <- renderPlot({
      re <- lme4::ranef(sleep_fit(), condVar = TRUE)$Subject
      pv <- attr(re, "postVar")
      se <- sqrt(t(apply(pv, 3, diag)))
      if (ncol(re) == 1) se <- matrix(se, ncol = 1)
      lab <- c(`(Intercept)` = "Baseline (ms)", Days = "Per day (ms)")
      df <- do.call(rbind, lapply(seq_len(ncol(re)), function(j)
        data.frame(subject = rownames(re), term = unname(lab[colnames(re)[j]]), est = re[, j], se = se[, j])))
      ord <- rownames(re)[order(re[, 1])]
      df$subject <- factor(df$subject, levels = ord)
      ggplot(df, aes(est, subject)) +
        geom_vline(xintercept = 0, colour = sv_col$ink2, linetype = "22") +
        geom_errorbarh(aes(xmin = est - 1.96 * se, xmax = est + 1.96 * se), height = 0, colour = sv_col$primary,
                       alpha = 0.6, linewidth = 0.8) +
        geom_point(colour = sv_col$primary, size = 2) +
        facet_wrap(~term, scales = "free_x") +
        labs(x = "Difference from average", y = NULL) +
        theme_statviz(10) + theme(strip.text = element_text(face = "bold", colour = sv_col$ink))
    })

    quiz_server(input, output, session, list(
      list(q = "Why can't we just run an ordinary regression on students from several schools?",
           choices = c("Regression needs at least 1000 rows",
                       "Students in the same school aren't independent, so the standard errors come out too small",
                       "Schools are categorical", "The slope would always be zero"),
           answer = 2,
           why = "Shared school effects make students in a school alike; treating them as independent overstates how much information you have."),
      list(q = "In a mixed model, 'partial pooling' means...",
           choices = c("throwing away small groups",
                       "fitting one line to everybody",
                       "each group gets its own estimate, pulled toward the overall average",
                       "fitting completely separate models"),
           answer = 3,
           why = "It's the compromise between complete pooling and no pooling."),
      list(q = "Which school's estimate gets shrunk the most toward the average?",
           choices = c("The largest school", "A school with 2 students", "The school closest to average",
                       "All schools are shrunk equally"),
           answer = 2,
           why = "Little data = noisy estimate = more weight on the overall average."),
      list(q = "An ICC of 0.30 means...",
           choices = c("30% of schools are significant", "30% of the variation in the outcome is between groups",
                       "the slope is 0.30", "30% of students were excluded"),
           answer = 2,
           why = "ICC = between-group variance / total variance."),
      list(q = "In y = (β0 + u0j) + (β1 + u1j)x + ε, the 'random effects' are...",
           choices = c("β0 and β1", "u0j and u1j", "ε", "x"),
           answer = 2,
           why = "The u's are group-specific deviations drawn from a distribution; the β's are the fixed (average) effects.")
    ))
  })
}
