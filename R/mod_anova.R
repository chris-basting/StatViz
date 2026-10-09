# ANOVA module -----------------------------------------------------------------
#
# Tabs:
#   1. Why not many t-tests?  - multiple comparisons
#   2. Between vs within      - sums of squares and the F ratio
#   3. The F distribution     - p-value from F
#   4. Which groups differ?   - Tukey post-hoc
#   5. Two factors            - two-way ANOVA and interaction
#   6. Real data
#   7. Quiz

anova_ui <- function(id) {
  ns <- NS(id)

  layout_sidebar(
    sidebar = sidebar(
      width = 300,
      show_on(ns, "why",
        sliderInput(ns("k_why"), "Number of groups", 2, 10, 5, step = 1),
        p(class = "small text-muted", "Each pair of groups is compared with its own t-test at α = 0.05.")
      ),
      show_on(ns, c("idea", "f", "posthoc"),
        h6("Simulated groups"),
        sliderInput(ns("k"), "Number of groups", 2, 6, 4, step = 1),
        sliderInput(ns("n"), "Observations per group", 3, 50, 10, step = 1),
        sliderInput(ns("spread"), "How different the true group means are", 0, 4, 1, step = 0.1),
        sliderInput(ns("sd"), "Noise within groups (SD)", 0.5, 6, 2, step = 0.1),
        actionButton(ns("resample"), "Draw a new sample", icon = icon("shuffle"),
                     class = "btn-outline-secondary btn-sm w-100")
      ),
      show_on(ns, "idea",
        hr(),
        radioButtons(ns("show_dev"), "Show distances", c(
          "Total: each point to the grand mean" = "total",
          "Between: group means to the grand mean" = "between",
          "Within: each point to its group mean" = "within"
        ), selected = "between")
      ),
      show_on(ns, "twoway",
        h6("Plant growth: fertilizer × light"),
        sliderInput(ns("eff_a"), "Effect of fertilizer", 0, 4, 2, step = 0.1),
        sliderInput(ns("eff_b"), "Effect of light", 0, 4, 1.5, step = 0.1),
        sliderInput(ns("inter"), "Interaction (fertilizer helps more in bright light)", -4, 4, 0, step = 0.1),
        sliderInput(ns("sd2"), "Noise (SD)", 0.5, 5, 1.5, step = 0.1),
        sliderInput(ns("n2"), "Plants per combination", 3, 30, 8, step = 1),
        actionButton(ns("resample2"), "New plants", icon = icon("shuffle"),
                     class = "btn-outline-secondary btn-sm w-100")
      ),
      show_on(ns, "real",
        selectInput(ns("dataset"), "Dataset",
                    setNames(names(anova_datasets), vapply(anova_datasets, `[[`, "", "name"))),
        uiOutput(ns("blurb")),
        uiOutput(ns("twoway_ui"))
      ),
      show_on(ns, "quiz", p("Test yourself on the ideas from the other tabs."))
    ),

    navset_underline(
      id = ns("tab"),

      # ==== 1. Why =============================================================
      nav_panel(
        "1 · Why not many t-tests?", value = "why",
        explain(withMathJax(
          p("With more than two groups, why not just run a t-test for every pair? Because every test has a",
            "5% chance of a false alarm, and those chances add up. With 5 groups there are 10 pairs,",
            "so the odds that at least one comparison looks 'significant' by pure chance are about 40%."),
          math_words(
            "$$P(\\text{at least one false alarm}) = 1 - (1 - 0.05)^{m}, \\qquad m = \\frac{k(k-1)}{2}$$",
            "the chance that all m comparisons stay quiet is 0.95 multiplied by itself m times; one minus that is",
            "the chance that at least one goes off. (This treats the tests as independent, so it's approximate.)",
            symbols = list("\\(k\\)" = "number of groups", "\\(m\\)" = "number of pairs of groups")
          ),
          p(strong("ANOVA"), " (analysis of variance) answers the overall question, 'are any of the group",
            "means different?', with a single test and a single 5% false-alarm rate.",
            "If the answer is yes, follow-up (post-hoc) tests that adjust for multiple comparisons say which ones.")
        )),
        layout_columns(
          col_widths = c(7, 5),
          card(card_header("False alarms pile up"),
               plotOutput(ns("fwer_plot"), height = "360px"),
               axes_note("number of groups.", "chance of at least one false alarm across all pairwise t-tests.",
                         "Dotted line: the 5% you'd expect from a single test. Orange dot: your number of groups.")),
          div(
            uiOutput(ns("fwer_boxes")),
            card(card_header("Check it by simulation"),
                 actionButton(ns("sim_fwer"), "Run 500 experiments with no real differences", icon = icon("play"),
                              class = "btn-outline-primary btn-sm"),
                 uiOutput(ns("fwer_sim")))
          )
        )
      ),

      # ==== 2. Between vs within ==============================================
      nav_panel(
        "2 · Between vs within", value = "idea",
        explain(withMathJax(
          p("ANOVA compares two kinds of spread. If the groups really differ, the group averages will be",
            "spread out ", em("more"), " than you'd expect from the scatter of individuals inside each group."),
          math_words(
            "$$\\underset{\\text{total}}{\\text{SS}_T} = \\underset{\\text{between groups}}{\\text{SS}_B} + \\underset{\\text{within groups}}{\\text{SS}_W}$$",
            "each point's squared distance from the overall average splits into two parts: how far its group's",
            "average is from the overall average (between), plus how far the point is from its own group's average",
            "(within). Use the sidebar to draw each kind of distance. It's the same 'total = explained + left over'",
            "idea as R² in regression."
          ),
          math_words(
            "$$F = \\frac{\\text{SS}_B / (k-1)}{\\text{SS}_W / (N-k)} = \\frac{\\text{MS}_B}{\\text{MS}_W}$$",
            "divide each sum of squares by its ", strong("degrees of freedom"), " to get a 'mean square' (an average",
            "spread), then take the ratio. If there are no real differences, both mean squares estimate the same",
            "noise and F is around 1. Real differences make the top bigger, so F grows.",
            symbols = list("\\(k\\)" = "number of groups", "\\(N\\)" = "total number of observations",
                           "\\(\\text{MS}\\)" = "mean square = SS ÷ degrees of freedom")
          )
        )),
        layout_columns(
          col_widths = c(7, 5),
          card(full_screen = TRUE, card_header("Groups, group means and the grand mean"),
               plotOutput(ns("dev_plot"), height = "430px"),
               axes_note("group.", "measured value.",
                         "Dots: observations. Black dashes: group averages. Dotted line: the grand (overall) average.",
                         "Orange segments: the distances selected in the sidebar.")),
          div(
            card(card_header("The ANOVA table"), uiOutput(ns("aov_table"))),
            uiOutput(ns("aov_boxes"))
          )
        )
      ),

      # ==== 3. F distribution =================================================
      nav_panel(
        "3 · The F distribution", value = "f",
        explain(
          p("When there are no real differences between groups, F still varies by chance, following the",
            strong("F distribution"), ". Its shape depends on the two degrees of freedom (number of groups and",
            "number of observations)."),
          p("Only ", em("large"), " values of F count as evidence of a difference (small F just means the group",
            "means are unusually close together), so the p-value is the area in the right-hand tail only."),
          p("Try: set 'How different the true group means are' to 0 and redraw a few samples. F hovers around 1",
            "and p jumps around between 0 and 1. That's what 'no effect' looks like.")
        ),
        layout_columns(
          col_widths = c(6, 6),
          card(card_header("Your F compared with the no-difference distribution"),
               plotOutput(ns("f_curve"), height = "340px"),
               axes_note("F statistic.", "how likely each F value is when all group means are equal.",
                         "Black line: your F. Orange area: p-value.")),
          card(card_header("500 experiments with no real differences"),
               actionButton(ns("sim_f"), "Run them", icon = icon("play"), class = "btn-outline-primary btn-sm mb-2",
                            width = "160px"),
               plotOutput(ns("f_hist"), height = "290px"),
               axes_note("F from each simulated experiment.", "number of experiments.",
                         "Orange: experiments with F at least as large as yours."))
        )
      ),

      # ==== 4. Post-hoc =======================================================
      nav_panel(
        "4 · Which groups differ?", value = "posthoc",
        explain(
          p("A significant F says ", em("some"), " group means differ, but not which. ",
            strong("Tukey's HSD"), " ('honestly significant difference') compares every pair while keeping the",
            "overall false-alarm rate at 5%. It does this by making each interval a bit wider than an ordinary",
            "t-test interval would be."),
          p("In the chart, a pair whose interval doesn't cross 0 is a clear difference. Compare the blue (Tukey)",
            "intervals with the gray (unadjusted) ones: the unadjusted ones are narrower, and that's exactly",
            "what lets false alarms slip through.")
        ),
        card(full_screen = TRUE, card_header("All pairwise differences"),
             plotOutput(ns("tukey_plot"), height = "440px"),
             axes_note("difference between the two group means (second minus first), with 95% interval.",
                       "pair of groups.",
                       "Blue: Tukey-adjusted intervals. Gray: unadjusted t intervals. Dashed line: no difference."))
      ),

      # ==== 5. Two-way ========================================================
      nav_panel(
        "5 · Two factors", value = "twoway",
        explain(withMathJax(
          p("Often two things vary at once, e.g. fertilizer (none / added) and light (low / medium / high).",
            strong("Two-way ANOVA"), " asks three questions at once:"),
          tags$ul(
            tags$li(strong("Main effect of fertilizer:"), " averaging over light levels, does fertilizer matter?"),
            tags$li(strong("Main effect of light:"), " averaging over fertilizer, does light matter?"),
            tags$li(strong("Interaction:"), " does the effect of fertilizer ", em("depend on"), " the light level?")
          ),
          p("In an ", strong("interaction plot"), " the lines are parallel when there's no interaction (fertilizer",
            "adds the same amount everywhere) and diverge or cross when there is one. Move the interaction slider",
            "and watch the lines fan out."),
          deeper("Show the model",
                 math_words("$$y = \\mu + \\alpha_{\\text{fert}} + \\beta_{\\text{light}} + (\\alpha\\beta)_{\\text{fert, light}} + \\varepsilon$$",
                            "each plant's growth = overall average + a bump for its fertilizer + a bump for its light",
                            "level + an extra bump for that particular combination (the interaction) + noise."))
        )),
        layout_columns(
          col_widths = c(7, 5),
          card(full_screen = TRUE, card_header("Interaction plot"),
               plotOutput(ns("inter_plot"), height = "400px"),
               axes_note("light level.", "plant growth.",
                         "Faint dots: individual plants. Lines join the average for each fertilizer at each light level.")),
          card(card_header("Two-way ANOVA table"), uiOutput(ns("two_table")))
        )
      ),

      # ==== 6. Real data ======================================================
      nav_panel(
        "6 · Real data", value = "real",
        explain(p(
          "Pick a dataset. One-factor datasets get a one-way ANOVA and Tukey comparisons. Two-factor datasets",
          "let you add the second factor and its interaction. ANOVA assumes roughly bell-shaped noise with similar",
          "spread in every group; if that's badly off (e.g. insect counts), see the Kruskal-Wallis test on the",
          "Wilcoxon page."
        )),
        layout_columns(
          col_widths = c(6, 6),
          card(card_header("The data"),
               plotOutput(ns("real_plot"), height = "380px"),
               uiOutput(ns("real_axes"))),
          card(card_header("ANOVA"), uiOutput(ns("real_table")))
        ),
        card(card_header("Tukey comparisons for the first factor"),
             plotOutput(ns("real_tukey"), height = "320px"),
             axes_note("difference between group means, with Tukey-adjusted 95% interval.", "pair of groups.",
                       "Intervals that don't cross the dashed zero line show a clear difference."))
      ),

      nav_panel("Quiz", value = "quiz", quiz_panel(ns))
    )
  )
}

anova_server <- function(id) {
  moduleServer(id, function(input, output, session) {

    # ---- 1. why -------------------------------------------------------------
    fwer <- function(k) 1 - 0.95^(k * (k - 1) / 2)

    output$fwer_plot <- renderPlot({
      df <- data.frame(k = 2:10, p = fwer(2:10))
      ggplot(df, aes(k, p)) +
        geom_hline(yintercept = 0.05, linetype = "dotted", colour = sv_col$ink2) +
        geom_line(colour = sv_col$primary, linewidth = 1.1) +
        geom_point(colour = sv_col$primary, size = 2.5) +
        annotate("point", x = input$k_why, y = fwer(input$k_why), colour = sv_col$accent, size = 5) +
        scale_x_continuous(breaks = 2:10) +
        scale_y_continuous(limits = c(0, 1), labels = function(v) paste0(100 * v, "%")) +
        labs(x = "Number of groups", y = "Chance of ≥ 1 false alarm") + theme_statviz()
    })

    output$fwer_boxes <- renderUI({
      k <- input$k_why
      layout_column_wrap(
        width = 1 / 2, fill = FALSE,
        sv_box("Pairs to compare", k * (k - 1) / 2, tone = "blue"),
        sv_box("Chance of ≥ 1 false alarm", sprintf("%.0f%%", 100 * fwer(k)), tone = "orange")
      )
    })

    fwer_res <- reactiveVal(NULL)
    observeEvent(input$k_why, fwer_res(NULL))
    observeEvent(input$sim_fwer, {
      k <- input$k_why
      hits <- vapply(seq_len(500), function(i) {
        d <- sim_k_groups(k, 10, 0, 1, 7000 + i)
        any(pairwise.t.test(d$y, d$group, p.adjust.method = "none")$p.value < 0.05, na.rm = TRUE)
      }, TRUE)
      fwer_res(mean(hits))
    })
    output$fwer_sim <- renderUI({
      r <- fwer_res()
      if (is.null(r)) return(p(class = "small text-muted mt-2",
                               "Each simulated experiment has groups with identical true means and runs every pairwise t-test."))
      p(class = "mt-2", strong(sprintf("%.0f%%", 100 * r)), " of the 500 experiments had at least one 'significant' pair,",
        "even though nothing was going on. (The formula slightly overstates it because the tests share data.)")
    })

    # ---- shared simulated groups (tabs 2-4) ---------------------------------
    seed <- reactiveVal(71)
    observeEvent(input$resample, seed(sample.int(1e6, 1)))
    dat <- reactive(sim_k_groups(input$k, input$n, input$spread, input$sd, seed()))
    fit <- reactive(aov(y ~ group, data = dat()))
    tab <- reactive({
      s <- summary(fit())[[1]]
      list(df = s$Df, ss = s$`Sum Sq`, ms = s$`Mean Sq`, F = s$`F value`[1], p = s$`Pr(>F)`[1])
    })

    output$dev_plot <- renderPlot({
      d <- dat(); gm <- mean(d$y); m <- tapply(d$y, d$group, mean)
      set.seed(4)
      d$xj <- as.integer(d$group) + runif(nrow(d), -0.25, 0.25)
      d$gmean <- m[d$group]
      g <- ggplot(d) + geom_hline(yintercept = gm, colour = sv_col$ink2, linetype = "dotted", linewidth = 0.8)
      g <- g + switch(
        input$show_dev,
        total = geom_segment(aes(xj, y, xend = xj, yend = gm), colour = sv_col$accent, alpha = 0.6),
        within = geom_segment(aes(xj, y, xend = xj, yend = gmean), colour = sv_col$accent, alpha = 0.6),
        between = geom_segment(data = data.frame(x = seq_along(m), m = m), aes(x, m, xend = x, yend = gm),
                               colour = sv_col$accent, linewidth = 2.5)
      )
      g + geom_point(aes(xj, y), colour = sv_col$primary, alpha = 0.7, size = 2.2) +
        geom_segment(data = data.frame(x = seq_along(m), m = m), aes(x = x - 0.35, xend = x + 0.35, y = m, yend = m),
                     colour = sv_col$ink, linewidth = 1.4) +
        annotate("text", x = 0.45, y = gm, label = "grand mean", hjust = 0, vjust = -0.5, size = 3.4, colour = sv_col$ink2) +
        scale_x_continuous(breaks = seq_along(m), labels = names(m)) +
        labs(x = "Group", y = "Value") + theme_statviz()
    })

    aov_table_ui <- function(rows, sentence) {
      tagList(
        tags$table(
          class = "table table-sm",
          tags$thead(tags$tr(tags$th("Source"), tags$th("df"), tags$th("SS"), tags$th("MS"), tags$th("F"), tags$th("p"))),
          tags$tbody(lapply(rows, function(r) tags$tr(lapply(r, tags$td))))
        ),
        sentence
      )
    }

    output$aov_table <- renderUI({
      t <- tab()
      aov_table_ui(
        list(list("Between groups", t$df[1], fmt(t$ss[1], 1), fmt(t$ms[1]), strong(fmt(t$F)), strong(fmt_p(t$p))),
             list("Within groups", t$df[2], fmt(t$ss[2], 1), fmt(t$ms[2]), "", ""),
             list(em("Total"), sum(t$df), fmt(sum(t$ss), 1), "", "", "")),
        tagList(
          tags$ul(class = "small mb-1",
            tags$li(strong("df"), " (degrees of freedom): groups − 1 for 'between', observations − groups for 'within'."),
            tags$li(strong("SS"), ": sum of squared distances, as drawn in the chart."),
            tags$li(strong("MS"), " = SS ÷ df, an average spread."),
            tags$li(strong("F"), " = MS between ÷ MS within.")),
          p(class = "small mb-0", strong("In words: "), p_sentence(t$p))
        )
      )
    })

    output$aov_boxes <- renderUI({
      t <- tab()
      layout_column_wrap(
        width = 1 / 2, fill = FALSE,
        sv_box("F", fmt(t$F), p("≈ 1 when there are no real differences"), tone = "blue"),
        sv_box("Share of variation between groups", sprintf("%.0f%%", 100 * t$ss[1] / sum(t$ss)),
               p("SS between ÷ SS total (called η², 'eta-squared')"), tone = "green")
      )
    })

    # ---- 3. F distribution --------------------------------------------------
    output$f_curve <- renderPlot({
      t <- tab()
      x <- seq(0.001, max(6, t$F * 1.2), length.out = 400)
      dens <- df(x, t$df[1], t$df[2])
      null_curve_plot(x, dens, t$F, "upper", xlab = "F", p = t$p, ymax = min(1.5, max(dens)))
    })

    f_sims <- reactiveVal(NULL)
    observeEvent(list(input$k, input$n, input$spread, input$sd, input$resample), f_sims(NULL), ignoreInit = TRUE)
    observeEvent(input$sim_f, {
      k <- input$k; n <- input$n; sd <- input$sd
      f_sims(vapply(seq_len(500), function(i) {
        d <- sim_k_groups(k, n, 0, sd, 9000 + i)
        summary(aov(y ~ group, data = d))[[1]]$`F value`[1]
      }, 0))
    })

    output$f_hist <- renderPlot({
      s <- f_sims()
      validate(need(!is.null(s), "Press 'Run them' to simulate."))
      f0 <- tab()$F
      df <- data.frame(F = s, extreme = s >= f0)
      ggplot(df, aes(F, fill = extreme)) +
        geom_histogram(binwidth = 0.2, boundary = 0, colour = sv_col$surface) +
        geom_vline(xintercept = f0, colour = sv_col$ink) +
        scale_fill_manual(values = c(`FALSE` = "#9ec5f4", `TRUE` = sv_col$accent), guide = "none") +
        annotate("label", x = Inf, y = Inf, hjust = 1.05, vjust = 1.3, size = 3.8, fill = sv_col$surface, label.size = 0,
                 label = sprintf("%d of 500 (%.1f%%) at least as large", sum(df$extreme), 100 * mean(df$extreme))) +
        coord_cartesian(xlim = c(0, max(6, f0 * 1.2))) +
        labs(x = "F", y = "Experiments") + theme_statviz()
    })

    # ---- 4. Tukey -----------------------------------------------------------
    tukey_df <- function(fit, d, y, g) {
      tk <- TukeyHSD(fit, g)[[1]]
      lv <- levels(d[[g]]); s <- sqrt(summary(fit)[[1]]$`Mean Sq`[nrow(summary(fit)[[1]])])
      dfres <- fit$df.residual
      out <- data.frame(pair = rownames(tk), diff = tk[, "diff"], lo = tk[, "lwr"], hi = tk[, "upr"], p = tk[, "p adj"])
      # unadjusted t intervals using the same pooled SD
      ns <- table(d[[g]])
      parts <- strsplit(out$pair, "-", fixed = TRUE)
      se <- vapply(parts, function(pp) s * sqrt(1 / ns[[pp[1]]] + 1 / ns[[pp[2]]]), 0)
      out$ulo <- out$diff - qt(0.975, dfres) * se
      out$uhi <- out$diff + qt(0.975, dfres) * se
      out$pair <- factor(out$pair, levels = rev(out$pair))
      out
    }

    tukey_plot <- function(df) {
      ggplot(df, aes(diff, pair)) +
        geom_vline(xintercept = 0, linetype = "22", colour = sv_col$ink2) +
        geom_errorbarh(aes(xmin = ulo, xmax = uhi), height = 0, colour = sv_col$muted, linewidth = 3, alpha = 0.5) +
        geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.25, colour = sv_col$primary, linewidth = 0.9) +
        geom_point(colour = sv_col$primary, size = 3) +
        labs(x = "Difference in means", y = NULL) + theme_statviz()
    }

    output$tukey_plot <- renderPlot(tukey_plot(tukey_df(fit(), dat(), "y", "group")))

    # ---- 5. two-way ---------------------------------------------------------
    seed2 <- reactiveVal(81)
    observeEvent(input$resample2, seed2(sample.int(1e6, 1)))
    dat2 <- reactive({
      set.seed(seed2())
      grid <- expand.grid(fert = c("None", "Added"), light = c("Low", "Medium", "High"))
      grid <- grid[rep(seq_len(nrow(grid)), each = input$n2), ]
      a <- as.integer(grid$fert == "Added"); b <- as.integer(grid$light) - 1   # 0, 1, 2
      grid$y <- 10 + input$eff_a * a + input$eff_b * b + input$inter * a * (b - 1) + rnorm(nrow(grid), sd = input$sd2)
      grid$fert <- factor(grid$fert, levels = c("None", "Added"))
      grid
    })

    output$inter_plot <- renderPlot({
      d <- dat2()
      m <- aggregate(y ~ fert + light, data = d, FUN = mean)
      cols <- c(None = sv_col$primary, Added = sv_col$accent)
      ggplot(d, aes(light, y, colour = fert)) +
        geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.3, seed = 1), alpha = 0.35, size = 1.8) +
        geom_line(data = m, aes(group = fert), linewidth = 1.3, position = position_dodge(0.3)) +
        geom_point(data = m, size = 4, position = position_dodge(0.3)) +
        scale_colour_manual(values = cols, name = "Fertilizer") +
        labs(x = "Light", y = "Growth") + theme_statviz()
    })

    two_rows <- function(fit, labels) {
      s <- summary(fit)[[1]]
      rn <- trimws(rownames(s))
      lapply(seq_len(nrow(s)), function(i) {
        lab <- if (rn[i] %in% names(labels)) labels[[rn[i]]] else rn[i]
        list(lab, s$Df[i], fmt(s$`Sum Sq`[i], 1), fmt(s$`Mean Sq`[i]),
             if (is.na(s$`F value`[i])) "" else strong(fmt(s$`F value`[i])),
             if (is.na(s$`Pr(>F)`[i])) "" else strong(fmt_p(s$`Pr(>F)`[i])))
      })
    }

    effect_sentences <- function(fit, labels) {
      s <- summary(fit)[[1]]; rn <- trimws(rownames(s))
      tags$ul(class = "small mb-0", lapply(which(!is.na(s$`Pr(>F)`)), function(i) {
        lab <- if (rn[i] %in% names(labels)) labels[[rn[i]]] else rn[i]
        tags$li(strong(lab), ": ", p_sentence(s$`Pr(>F)`[i]))
      }))
    }

    output$two_table <- renderUI({
      f <- aov(y ~ fert * light, data = dat2())
      labs <- list(fert = "Fertilizer", light = "Light", `fert:light` = "Fertilizer × light (interaction)",
                   Residuals = "Within (residual)")
      aov_table_ui(two_rows(f, labs), effect_sentences(f, labs))
    })

    # ---- 6. real data -------------------------------------------------------
    ds <- reactive(anova_datasets[[input$dataset]])
    output$blurb <- renderUI(p(class = "small text-muted", ds()$blurb))
    output$twoway_ui <- renderUI({
      if (!is.null(ds()$b)) tagList(
        checkboxInput(session$ns("use_b"), sprintf("Include second factor (%s)", ds()$b), TRUE),
        checkboxInput(session$ns("use_int"), "Include interaction", TRUE)
      )
    })
    use_b <- reactive(!is.null(ds()$b) && isTRUE(input$use_b))

    real_fit <- reactive({
      d <- ds()$data
      rhs <- if (use_b()) paste(ds()$a, if (isTRUE(input$use_int)) "*" else "+", ds()$b) else ds()$a
      aov(as.formula(paste(ds()$y, "~", rhs)), data = d)
    })

    output$real_plot <- renderPlot({
      d <- ds()$data
      if (use_b()) {
        m <- aggregate(as.formula(paste(ds()$y, "~", ds()$a, "+", ds()$b)), data = d, FUN = mean)
        bl <- levels(factor(d[[ds()$b]]))
        cols <- setNames(sv_groups[seq_along(bl)], bl)
        ggplot(d, aes(.data[[ds()$a]], .data[[ds()$y]], colour = .data[[ds()$b]])) +
          geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.3, seed = 1), alpha = 0.35) +
          geom_line(data = m, aes(group = .data[[ds()$b]]), linewidth = 1.2, position = position_dodge(0.3)) +
          geom_point(data = m, size = 3.5, position = position_dodge(0.3)) +
          scale_colour_manual(values = cols) + theme_statviz()
      } else {
        d$g <- factor(d[[ds()$a]])
        group_points_plot(d, y = ds()$y, group = "g", ylab = ds()$y) + labs(x = ds()$a)
      }
    })

    output$real_axes <- renderUI({
      if (use_b()) axes_note(paste0(ds()$a, "."), paste0(ds()$y, "."),
                             sprintf("Faint dots: observations. Lines join the average for each %s.", ds()$b))
      else axes_note(paste0(ds()$a, " (the groups)."), paste0(ds()$y, "."),
                     "Dots: observations. Black: group average with 95% CI.")
    })

    output$real_table <- renderUI({
      f <- real_fit()
      labs <- list(Residuals = "Within (residual)")
      if (use_b()) labs[[paste0(ds()$a, ":", ds()$b)]] <- sprintf("%s × %s (interaction)", ds()$a, ds()$b)
      aov_table_ui(two_rows(f, labs), effect_sentences(f, labs))
    })

    output$real_tukey <- renderPlot({
      d <- ds()$data; d[[ds()$a]] <- factor(d[[ds()$a]])
      tukey_plot(tukey_df(real_fit(), d, ds()$y, ds()$a))
    })

    quiz_server(input, output, session, list(
      list(q = "Why use ANOVA instead of t-tests on every pair of groups?",
           choices = c("t-tests can't handle more than 10 observations",
                       "Many separate tests inflate the chance of at least one false alarm",
                       "ANOVA tells you exactly which groups differ", "ANOVA doesn't need any assumptions"),
           answer = 2, why = "With m tests at 5% each, the chance of some false alarm grows quickly. ANOVA asks one overall question."),
      list(q = "F is the ratio of...",
           choices = c("the largest mean to the smallest mean",
                       "spread between group means to spread within groups",
                       "the number of groups to the number of observations", "two p-values"),
           answer = 2, why = "F = MS between ÷ MS within. Real differences inflate the between part."),
      list(q = "If all true group means are equal, F will usually be about...",
           choices = c("0", "1", "the number of groups", "100"),
           answer = 2, why = "Both mean squares then estimate the same noise, so their ratio hovers around 1."),
      list(q = "ANOVA gives p < 0.05. What do you know?",
           choices = c("Every group differs from every other group", "At least one group mean differs from the others",
                       "The first group is the largest", "The data are normally distributed"),
           answer = 2, why = "ANOVA is an overall test; use post-hoc comparisons like Tukey's HSD to find which."),
      list(q = "In an interaction plot, the lines are clearly not parallel. This suggests...",
           choices = c("no effect of either factor", "the effect of one factor depends on the level of the other",
                       "the data are skewed", "the sample is too small"),
           answer = 2, why = "Non-parallel lines mean the gap between levels of one factor changes across the other: an interaction.")
    ))
  })
}
