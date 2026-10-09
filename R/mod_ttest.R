# t-test module ----------------------------------------------------------------
#
# Tabs:
#   1. Signal vs noise     - the t statistic
#   2. Where p-values come from - null distribution, simulation
#   3. Paired vs unpaired  - why pairing helps; one-sample idea
#   4. Errors & power      - alpha, beta, power
#   5. Real data
#   6. Quiz

ttest_ui <- function(id) {
  ns <- NS(id)

  layout_sidebar(
    sidebar = sidebar(
      width = 300,
      show_on(ns, c("signal", "null"),
        h6("Simulated experiment"),
        p(class = "small text-muted",
          "Two groups, e.g. plants with and without a fertilizer. You choose the true difference."),
        sliderInput(ns("n"), "Observations per group", 3, 100, 12, step = 1),
        sliderInput(ns("diff"), "True difference between groups", 0, 4, 1, step = 0.1),
        sliderInput(ns("sd"), "Noise (SD within each group)", 0.5, 5, 1.5, step = 0.1),
        actionButton(ns("resample"), "Run the experiment again", icon = icon("shuffle"),
                     class = "btn-outline-secondary btn-sm w-100")
      ),
      show_on(ns, "paired",
        h6("Before / after on the same people"),
        sliderInput(ns("np"), "Number of people", 3, 50, 10, step = 1),
        sliderInput(ns("change"), "True change (after − before)", -5, 5, 2, step = 0.5),
        sliderInput(ns("sd_person"), "How much people differ from each other (SD)", 0, 20, 10, step = 1),
        sliderInput(ns("sd_noise"), "Measurement noise (SD)", 0.5, 10, 3, step = 0.5),
        actionButton(ns("resample_p"), "New people", icon = icon("shuffle"),
                     class = "btn-outline-secondary btn-sm w-100")
      ),
      show_on(ns, "power",
        sliderInput(ns("d"), "True effect size (difference ÷ SD)", 0, 2, 0.5, step = 0.05),
        sliderInput(ns("n_pow"), "Observations per group", 2, 200, 20, step = 1),
        radioButtons(ns("alpha"), "Significance level (α)", c("0.10" = 0.10, "0.05" = 0.05, "0.01" = 0.01),
                     selected = 0.05, inline = TRUE)
      ),
      show_on(ns, "real",
        selectInput(ns("dataset"), "Dataset",
                    setNames(names(tt_datasets), vapply(tt_datasets, `[[`, "", "name"))),
        uiOutput(ns("blurb")),
        selectizeInput(ns("levels"), "Compare these two groups", NULL, multiple = TRUE,
                       options = list(maxItems = 2)),
        uiOutput(ns("paired_ui")),
        checkboxInput(ns("welch"), "Don't assume equal variances (Welch)", TRUE)
      ),
      show_on(ns, "quiz", p("Test yourself on the ideas from the other tabs."))
    ),

    navset_underline(
      id = ns("tab"),

      # ==== 1. Signal vs noise =================================================
      nav_panel(
        "1 · Signal vs noise", value = "signal",
        explain(withMathJax(
          p(strong("The question:"), "the treatment group's average is higher than the control's. Is that a",
            "real effect, or could it just be luck in who ended up in which group?"),
          p("The answer depends on two things: how big the difference is (", strong("signal"), ") and how much",
            "the averages would bounce around by chance (", strong("noise"), "). The ", strong("t statistic"),
            " is simply their ratio."),
          math_words(
            "$$t = \\frac{\\bar{y}_T - \\bar{y}_C}{\\text{SE}} \\qquad \\text{SE} = \\sqrt{\\frac{s_T^2}{n_T} + \\frac{s_C^2}{n_C}}$$",
            "difference between the two group averages, divided by its ", strong("standard error"),
            " (how much that difference would typically wobble if we re-ran the experiment).",
            "The SE is small when individuals vary little (small s) or when the groups are big (large n).",
            "A t around 0 means 'no visible signal'; beyond about ±2 the signal stands out from the noise.",
            symbols = list("\\(\\bar{y}_T, \\bar{y}_C\\)" = "group averages",
                           "\\(s_T, s_C\\)" = "group standard deviations",
                           "\\(n_T, n_C\\)" = "group sizes")
          ),
          p("Try it: keep the true difference fixed and raise the noise, or lower the sample size.",
            "The difference in averages stays about the same, but t shrinks.")
        )),
        layout_columns(
          col_widths = c(7, 5),
          card(full_screen = TRUE, card_header("The two groups"),
               plotOutput(ns("signal_plot"), height = "420px"),
               axes_note("group.", "measured value.",
                         "Dots: individual observations. Black diamond and bar: group average with its 95% confidence interval.")),
          div(uiOutput(ns("signal_boxes")))
        )
      ),

      # ==== 2. p-values ========================================================
      nav_panel(
        "2 · Where p-values come from", value = "null",
        explain(
          p("To judge whether our t is 'big', we compare it with what t looks like when there is ",
            strong("no real difference"), " (the ", strong("null hypothesis"), "). Even then, t won't be",
            "exactly 0, because of chance. The bell-shaped ", strong("t distribution"), " describes how far",
            "it typically wanders."),
          p("The ", strong("p-value"), " is the shaded area: the chance of a t at least as far from 0 as ours",
            em("if there were truly no difference"), ". Small p = our result would be surprising under 'no",
            "difference', which counts as evidence against it. The usual cut-off is 0.05."),
          p("Press the button to run 1,000 experiments where the true difference really is 0, and see how",
            "often chance alone produces a t as extreme as yours."),
          div(class = "small",
              strong("A p-value is NOT:"),
              tags$ul(class = "mb-0",
                tags$li("the probability that the null hypothesis is true;"),
                tags$li("the probability your result was 'due to chance';"),
                tags$li("a measure of how big or important the effect is (a tiny effect can have a tiny p with enough data).")))
        ),
        layout_columns(
          col_widths = c(6, 6),
          card(card_header("Your t compared with the null distribution"),
               plotOutput(ns("null_curve"), height = "340px"),
               axes_note("t statistic.", "how likely each t value is when there's no real difference.",
                         "Black line: your t (dashed: the same distance on the other side).",
                         "Orange area = p-value.")),
          card(card_header("1,000 experiments with no real difference"),
               actionButton(ns("sim_null"), "Run 1,000 'no difference' experiments", icon = icon("play"),
                            class = "btn-outline-primary btn-sm mb-2", width = "320px"),
               plotOutput(ns("null_hist"), height = "290px"),
               axes_note("t statistic from each simulated experiment.", "number of experiments.",
                         "Orange bars: experiments at least as extreme as yours. Their share should match the p-value."))
        )
      ),

      # ==== 3. Paired ==========================================================
      nav_panel(
        "3 · Paired vs unpaired", value = "paired",
        explain(withMathJax(
          p("If the same people are measured twice (before and after), the two sets of numbers are",
            strong("paired"), ". People differ a lot from each other, but each person's own change can be",
            "quite consistent. A ", strong("paired t-test"), " works on each person's change, so the",
            "person-to-person differences cancel out."),
          math_words(
            "$$t = \\frac{\\bar{d}}{s_d / \\sqrt{n}}$$",
            "compute each person's difference (after − before), then ask whether the average difference is far",
            "from 0 compared with how much the differences vary. This is the same as a ",
            strong("one-sample t-test"), " on the differences.",
            symbols = list("\\(d\\)" = "one person's after − before",
                           "\\(\\bar{d}, s_d\\)" = "average and SD of the differences",
                           "\\(n\\)" = "number of people")
          ),
          p("Try: raise 'How much people differ from each other'. The unpaired test loses its power",
            "(the groups overlap), but the paired test doesn't care.")
        )),
        layout_columns(
          col_widths = c(6, 6),
          card(full_screen = TRUE, card_header("Each person, before and after"),
               plotOutput(ns("paired_lines"), height = "380px"),
               axes_note("time point.", "measurement.",
                         "Each gray line joins one person's two measurements. Blue: lines going up; orange: going down.")),
          card(full_screen = TRUE, card_header("Each person's change"),
               plotOutput(ns("paired_diffs"), height = "380px"),
               axes_note("after − before.", "",
                         "Each dot is one person's change. Black: average change with 95% CI. Dashed line: no change."))
        ),
        uiOutput(ns("paired_table"))
      ),

      # ==== 4. Power ===========================================================
      nav_panel(
        "4 · Errors & power", value = "power",
        explain(
          p("Any test can be wrong in two ways:"),
          tags$ul(
            tags$li(strong("Type I error (false alarm):"), " concluding there's a difference when there isn't.",
                    "Its probability is the ", strong("significance level α"), " you pick (usually 0.05)."),
            tags$li(strong("Type II error (miss):"), " failing to detect a difference that really exists.")
          ),
          p(strong("Power"), " is the chance of detecting a real effect: 1 − P(miss). It grows with the",
            "size of the effect and with the sample size, and shrinks if you demand a stricter α.",
            "Studies are usually planned to have at least 80% power."),
          p("The effect size here is the difference divided by the SD (", strong("Cohen's d"), "):",
            "0.2 is often called small, 0.5 medium, 0.8 large.")
        ),
        layout_columns(
          col_widths = c(7, 5),
          card(full_screen = TRUE, card_header("What t looks like with and without a real effect"),
               plotOutput(ns("power_curves"), height = "380px"),
               axes_note("t statistic.", "how likely each t value is.",
                         "Gray: t when there's no effect. Blue: t when the effect is real. Dashed lines: the cut-offs",
                         "for 'significant'. Orange = false-alarm rate (α). Blue shading = power.")),
          div(
            uiOutput(ns("power_boxes")),
            card(card_header("Power vs sample size"),
                 plotOutput(ns("power_n"), height = "230px"),
                 axes_note("observations per group.", "power.",
                           "Dotted line: the common 80% target. Orange dot: your current choice."))
          )
        )
      ),

      # ==== 5. Real data =======================================================
      nav_panel(
        "5 · Real data", value = "real",
        explain(p(
          "Pick a dataset and two groups. The ", strong("Welch"), " version (default) doesn't assume the two",
          "groups have the same spread, which is the safer choice. The classic 'Student' version does",
          "assume equal spreads."
        )),
        layout_columns(
          col_widths = c(6, 6),
          card(card_header("The two groups"),
               plotOutput(ns("real_plot"), height = "380px"),
               axes_note("group.", "the measured value.",
                         "Dots: observations (gray lines join paired measurements). Black: average and 95% CI.")),
          card(card_header("Result"), uiOutput(ns("real_result")))
        )
      ),

      nav_panel("Quiz", value = "quiz", quiz_panel(ns))
    )
  )
}

ttest_server <- function(id) {
  moduleServer(id, function(input, output, session) {

    seed <- reactiveVal(51)
    observeEvent(input$resample, seed(sample.int(1e6, 1)))
    dat <- reactive(sim_two_groups(input$n, input$diff, input$sd, seed()))
    tt <- reactive(t.test(y ~ group, data = dat()))   # Welch; estimate order is Control, Treatment
    tval <- reactive(-unname(tt()$statistic))           # Treatment − Control

    # ---- 1. signal ----------------------------------------------------------
    output$signal_plot <- renderPlot({
      group_points_plot(dat(), ylab = "Measured value", colours = two_cols)
    })

    output$signal_boxes <- renderUI({
      d <- dat(); m <- tapply(d$y, d$group, mean); s <- tapply(d$y, d$group, sd); n <- tapply(d$y, d$group, length)
      se <- sqrt(sum(s^2 / n))
      tagList(
        sv_box("Signal: difference in averages", fmt(m[2] - m[1]),
               p(sprintf("Treatment %s − Control %s (true difference %s)", fmt(m[2]), fmt(m[1]), fmt(input$diff))),
               tone = "blue"),
        sv_box("Noise: standard error", fmt(se),
               p("typical chance wobble in that difference"), tone = "gray"),
        sv_box("t = signal ÷ noise", fmt(tval()),
               p(p_lab(tt()$p.value)), tone = if (tt()$p.value < 0.05) "green" else "orange")
      )
    })

    # ---- 2. null distribution -----------------------------------------------
    output$null_curve <- renderPlot({
      df <- tt()$parameter; t0 <- tval()
      lim <- max(4.5, abs(t0) + 1)
      x <- seq(-lim, lim, length.out = 400)
      null_curve_plot(x, dt(x, df), t0, "two", xlab = "t", p = tt()$p.value)
    })

    null_sims <- reactiveVal(NULL)
    observeEvent(list(input$n, input$diff, input$sd, input$resample), null_sims(NULL), ignoreInit = TRUE)
    observeEvent(input$sim_null, {
      n <- input$n; sd <- input$sd
      null_sims(vapply(seq_len(1000), function(i) {
        a <- rnorm(n, 10, sd); b <- rnorm(n, 10, sd)
        (mean(b) - mean(a)) / sqrt(var(a) / n + var(b) / n)
      }, 0))
    })

    output$null_hist <- renderPlot({
      s <- null_sims()
      validate(need(!is.null(s), "Press the button to run the simulation."))
      t0 <- abs(tval())
      df <- data.frame(t = s, extreme = abs(s) >= t0)
      lim <- max(4.5, t0 + 1)
      ggplot(df, aes(t, fill = extreme)) +
        geom_histogram(binwidth = 0.25, boundary = 0, colour = sv_col$surface) +
        geom_vline(xintercept = c(-t0, t0), colour = sv_col$ink, linetype = c("22", "solid")) +
        scale_fill_manual(values = c(`FALSE` = "#9ec5f4", `TRUE` = sv_col$accent), guide = "none") +
        annotate("label", x = -lim * 0.95, y = Inf, vjust = 1.3, hjust = 0, size = 3.8,
                 label = sprintf("%d of 1,000 (%.1f%%) as extreme as yours", sum(df$extreme), 100 * mean(df$extreme)),
                 fill = sv_col$surface, label.size = 0, colour = sv_col$ink) +
        coord_cartesian(xlim = c(-lim, lim)) +
        labs(x = "t", y = "Experiments") + theme_statviz()
    })

    # ---- 3. paired ----------------------------------------------------------
    seed_p <- reactiveVal(61)
    observeEvent(input$resample_p, seed_p(sample.int(1e6, 1)))
    pdat <- reactive(sim_paired(input$np, input$change, input$sd_person, input$sd_noise, seed_p()))

    output$paired_lines <- renderPlot({
      d <- pdat()
      long <- data.frame(id = rep(d$id, 2), time = factor(rep(c("Before", "After"), each = nrow(d)),
                                                          levels = c("Before", "After")),
                         y = c(d$before, d$after), up = rep(d$after > d$before, 2))
      ggplot(long, aes(time, y, group = id)) +
        geom_line(aes(colour = up), linewidth = 0.8, alpha = 0.8) +
        geom_point(colour = sv_col$ink2, size = 2) +
        scale_colour_manual(values = c(`TRUE` = sv_col$primary, `FALSE` = sv_col$accent), guide = "none") +
        labs(x = NULL, y = "Measurement") + theme_statviz()
    })

    output$paired_diffs <- renderPlot({
      d <- pdat(); dd <- data.frame(g = "Change", y = d$after - d$before)
      group_points_plot(dd, group = "g", ylab = "After − before") +
        geom_hline(yintercept = 0, colour = sv_col$ink2, linetype = "22") +
        coord_flip() + labs(x = NULL)
    })

    output$paired_table <- renderUI({
      d <- pdat()
      pt <- t.test(d$after, d$before, paired = TRUE)
      ut <- t.test(d$after, d$before)
      row <- function(name, t) tags$tr(tags$td(name), tags$td(fmt(t$estimate[1] - if (length(t$estimate) > 1) t$estimate[2] else 0)),
                                     tags$td(fmt(t$stderr)), tags$td(fmt(t$statistic)), tags$td(strong(fmt_p(t$p.value))))
      card(
        card_header("Same data, two tests"),
        tags$table(
          class = "table table-sm",
          tags$thead(tags$tr(tags$th("Test"), tags$th("Estimated change"), tags$th("SE"), tags$th("t"), tags$th("p-value"))),
          tags$tbody(row("Unpaired (treats before and after as unrelated groups)", ut),
                     row("Paired (uses each person's own change)", pt))
        ),
        p(class = "small mb-0",
          "Both tests estimate the same change, but the paired test's standard error is",
          if (pt$stderr < ut$stderr) sprintf(" %.1f times smaller,", ut$stderr / pt$stderr) else " not smaller here,",
          "because comparing each person with themselves removes the person-to-person variation.")
      )
    })

    # ---- 4. power -----------------------------------------------------------
    pw <- reactive({
      n <- input$n_pow; d <- input$d; a <- as.numeric(input$alpha)
      df <- 2 * n - 2; ncp <- d * sqrt(n / 2); crit <- qt(1 - a / 2, df)
      power <- if (d == 0) a else 1 - pt(crit, df, ncp) + pt(-crit, df, ncp)
      list(n = n, d = d, a = a, df = df, ncp = ncp, crit = crit, power = power)
    })

    output$power_curves <- renderPlot({
      w <- pw()
      x <- seq(-4.5, max(4.5, w$ncp + 4), length.out = 500)
      null <- dt(x, w$df); alt <- dt(x, w$df, w$ncp)
      df <- data.frame(x, null, alt)
      ggplot(df, aes(x)) +
        geom_area(aes(y = ifelse(abs(x) >= w$crit, alt, 0)), fill = sv_col$primary, alpha = 0.3) +
        geom_area(aes(y = ifelse(abs(x) >= w$crit, null, 0)), fill = sv_col$accent, alpha = 0.6) +
        geom_line(aes(y = null), colour = sv_col$muted, linewidth = 1) +
        geom_line(aes(y = alt), colour = sv_col$primary, linewidth = 1.2) +
        geom_vline(xintercept = c(-w$crit, w$crit), linetype = "22", colour = sv_col$ink) +
        annotate("text", x = 0, y = max(null) * 1.04, label = "no effect", colour = sv_col$ink2, vjust = 0) +
        annotate("text", x = w$ncp, y = max(alt) * 1.04 + 0.02, label = "real effect", colour = sv_col$primary,
                 vjust = 0, fontface = "bold") +
        scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
        labs(x = "t", y = "Probability density") + theme_statviz()
    })

    output$power_boxes <- renderUI({
      w <- pw()
      layout_column_wrap(
        width = 1 / 3, fill = FALSE,
        sv_box("Power", sprintf("%.0f%%", 100 * w$power), p("chance of detecting the effect"),
               tone = if (w$power >= 0.8) "green" else "blue"),
        sv_box("Miss rate", sprintf("%.0f%%", 100 * (1 - w$power)), p("Type II error"), tone = "gray"),
        sv_box("False alarm rate", sprintf("%.0f%%", 100 * w$a), p("Type I error (α)"), tone = "orange")
      )
    })

    output$power_n <- renderPlot({
      w <- pw()
      ns <- 2:200
      pw_n <- vapply(ns, function(n) {
        df <- 2 * n - 2; crit <- qt(1 - w$a / 2, df); ncp <- w$d * sqrt(n / 2)
        1 - pt(crit, df, ncp) + pt(-crit, df, ncp)
      }, 0)
      ggplot(data.frame(n = ns, p = pw_n), aes(n, p)) +
        geom_hline(yintercept = 0.8, linetype = "dotted", colour = sv_col$ink2) +
        geom_line(colour = sv_col$primary, linewidth = 1.1) +
        annotate("point", x = w$n, y = w$power, colour = sv_col$accent, size = 4) +
        scale_y_continuous(limits = c(0, 1), labels = function(v) paste0(100 * v, "%")) +
        labs(x = "Observations per group", y = "Power") + theme_statviz(11)
    })

    # ---- 5. real data -------------------------------------------------------
    ds <- reactive(tt_datasets[[input$dataset]])
    output$blurb <- renderUI(p(class = "small text-muted", ds()$blurb))
    observeEvent(ds(), {
      lv <- levels(factor(ds()$data[[ds()$group]]))
      updateSelectizeInput(session, "levels", choices = lv, selected = lv[1:2])
    })
    output$paired_ui <- renderUI({
      if (!is.null(ds()$id)) checkboxInput(session$ns("paired"), "Paired (the same subjects in both groups)", TRUE)
    })

    real_sub <- reactive({
      d <- ds()$data; lv <- input$levels
      validate(need(length(lv) == 2, "Pick exactly two groups."))
      d <- d[d[[ds()$group]] %in% lv, ]
      d$g <- factor(d[[ds()$group]], levels = lv)
      d
    })
    is_paired <- reactive(!is.null(ds()$id) && isTRUE(input$paired))

    real_test <- reactive({
      d <- real_sub(); y <- ds()$y; lv <- levels(d$g)
      if (is_paired()) {
        d <- d[order(d[[ds()$id]]), ]
        t.test(d[[y]][d$g == lv[2]], d[[y]][d$g == lv[1]], paired = TRUE)
      } else {
        t.test(d[[y]][d$g == lv[2]], d[[y]][d$g == lv[1]], var.equal = !isTRUE(input$welch))
      }
    })

    output$real_plot <- renderPlot({
      d <- real_sub(); y <- ds()$y
      g <- group_points_plot(d, y = y, group = "g", ylab = y, jitter = if (is_paired()) 0 else 0.12)
      if (is_paired()) g <- g + geom_line(aes(group = .data[[ds()$id]]), colour = sv_col$muted, alpha = 0.6)
      g
    })

    output$real_result <- renderUI({
      t <- real_test(); lv <- input$levels; y <- ds()$y
      est <- if (is_paired()) t$estimate else t$estimate[1] - t$estimate[2]
      kind <- if (is_paired()) "Paired t-test" else if (isTRUE(input$welch)) "Welch two-sample t-test" else "Student two-sample t-test"
      tagList(
        p(class = "text-muted small", kind),
        tags$table(
          class = "table table-sm",
          tags$tbody(
            tags$tr(tags$td(sprintf("Difference (%s − %s)", lv[2], lv[1])), tags$td(strong(fmt(est)))),
            tags$tr(tags$td("95% confidence interval"), tags$td(sprintf("%s to %s", fmt(t$conf.int[1]), fmt(t$conf.int[2])))),
            tags$tr(tags$td("t"), tags$td(fmt(t$statistic))),
            tags$tr(tags$td("Degrees of freedom"), tags$td(fmt(t$parameter, 1))),
            tags$tr(tags$td("p-value"), tags$td(strong(fmt_p(t$p.value))))
          )
        ),
        p(sprintf("In words: on average, %s is %s %s than %s", lv[2], fmt(abs(est)),
                  if (est >= 0) "higher" else "lower", lv[1]),
          sprintf("(%s), and we're fairly confident the true difference is between %s and %s.",
                  y, fmt(t$conf.int[1]), fmt(t$conf.int[2])),
          p_sentence(t$p.value))
      )
    })

    quiz_server(input, output, session, list(
      list(q = "The t statistic is...",
           choices = c("the difference in averages divided by its standard error",
                       "the probability the groups are equal", "the sum of the two averages",
                       "the standard deviation of the data"),
           answer = 1, why = "t = signal ÷ noise: how many standard errors apart the two averages are."),
      list(q = "Which change makes t larger, if the true difference stays the same?",
           choices = c("More noise within groups", "Fewer observations", "More observations", "None: t only depends on the difference"),
           answer = 3, why = "More data shrinks the standard error, so the same difference stands out more."),
      list(q = "p = 0.03 means...",
           choices = c("there's a 3% chance the null hypothesis is true",
                       "if there were no real difference, a result this extreme would happen about 3% of the time",
                       "the effect is small", "97% of the data support the treatment"),
           answer = 2, why = "A p-value is computed assuming no difference; it says how surprising your data would be under that assumption."),
      list(q = "When should you use a paired t-test?",
           choices = c("When the groups have equal size", "When the same subjects are measured in both conditions",
                       "When the data are skewed", "Always"),
           answer = 2, why = "Pairing compares each subject with themselves, removing person-to-person variation."),
      list(q = "Power is...",
           choices = c("the chance of a false alarm", "the chance of detecting a real effect",
                       "the size of the effect", "the p-value"),
           answer = 2, why = "Power = 1 − P(Type II error). It grows with effect size and sample size.")
    ))
  })
}
