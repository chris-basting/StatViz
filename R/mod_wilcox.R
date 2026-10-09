# Wilcoxon / rank-based tests module ---------------------------------------------
#
# Tabs:
#   1. Ranks instead of values - rank-sum statistic; effect of an outlier
#   2. Shuffling labels        - exact null distribution and permutation
#   3. t-test vs Wilcoxon      - false alarms and power under different shapes
#   4. Paired: signed-rank
#   5. More groups: Kruskal-Wallis
#   6. Real data
#   7. Quiz

wilcox_ui <- function(id) {
  ns <- NS(id)

  layout_sidebar(
    sidebar = sidebar(
      width = 300,
      show_on(ns, c("ranks", "null"),
        h6("Two small groups"),
        sliderInput(ns("n"), "Observations per group", 3, 15, 7, step = 1),
        sliderInput(ns("shift"), "True difference", 0, 4, 1.5, step = 0.1),
        selectInput(ns("shape"), "Shape of the noise", shape_choices),
        sliderInput(ns("outlier"), "Push one Treatment value far out by", 0, 40, 0, step = 1),
        actionButton(ns("resample"), "Draw a new sample", icon = icon("shuffle"),
                     class = "btn-outline-secondary btn-sm w-100")
      ),
      show_on(ns, "compare",
        selectInput(ns("shape_c"), "Shape of the noise", shape_choices, selected = "heavy"),
        sliderInput(ns("n_c"), "Observations per group", 5, 60, 20, step = 1),
        sliderInput(ns("shift_c"), "True difference (for the power test)", 0.2, 3, 1, step = 0.1),
        actionButton(ns("run_c"), "Run 500 experiments", icon = icon("play"), class = "btn-primary w-100")
      ),
      show_on(ns, "signed",
        sliderInput(ns("np"), "Number of people", 5, 30, 12, step = 1),
        sliderInput(ns("change"), "True average change", -3, 3, 1, step = 0.1),
        sliderInput(ns("sd_p"), "Noise in the changes (SD)", 0.5, 5, 2, step = 0.1),
        selectInput(ns("shape_p"), "Shape of the noise", shape_choices),
        actionButton(ns("resample_p"), "New people", icon = icon("shuffle"),
                     class = "btn-outline-secondary btn-sm w-100")
      ),
      show_on(ns, "kw",
        sliderInput(ns("k"), "Number of groups", 3, 5, 3, step = 1),
        sliderInput(ns("n_k"), "Observations per group", 4, 20, 8, step = 1),
        sliderInput(ns("spread"), "How different the groups are", 0, 3, 1, step = 0.1),
        selectInput(ns("shape_k"), "Shape of the noise", shape_choices, selected = "skewed"),
        actionButton(ns("resample_k"), "Draw a new sample", icon = icon("shuffle"),
                     class = "btn-outline-secondary btn-sm w-100")
      ),
      show_on(ns, "real",
        selectInput(ns("dataset"), "Dataset",
                    setNames(names(wx_datasets), vapply(wx_datasets, `[[`, "", "name"))),
        uiOutput(ns("blurb"))
      ),
      show_on(ns, "quiz", p("Test yourself on the ideas from the other tabs."))
    ),

    navset_underline(
      id = ns("tab"),

      # ==== 1. Ranks ===========================================================
      nav_panel(
        "1 · Ranks instead of values", value = "ranks",
        explain(withMathJax(
          p("The t-test compares ", em("averages"), ", so a single extreme value can drag an average (and inflate",
            "the standard deviation) a long way. The ", strong("Wilcoxon rank-sum test"), " (also called the",
            strong("Mann-Whitney U test"), ") sidesteps this by throwing away the actual values and keeping only",
            "their ", strong("order"), "."),
          p("Line up all observations from both groups, smallest to largest, and number them 1, 2, 3, ...",
            "These are the ", strong("ranks"), ". If the Treatment group tends to have larger values, its",
            "observations will collect the high ranks."),
          math_words(
            "$$W = \\sum (\\text{ranks of the Treatment group}) \\qquad \\text{expected if no difference} = \\frac{n_T (N+1)}{2}$$",
            "add up the Treatment group's ranks. If the groups were really the same, Treatment would get a",
            "random half of the ranks and its total would land near the expected value. A total far from that is",
            "evidence of a difference.",
            symbols = list("\\(n_T\\)" = "size of the Treatment group", "\\(N\\)" = "total number of observations")
          ),
          p("Try: push one Treatment value far out with the sidebar slider. The t-test's p-value gets",
            em("worse"), " (the outlier inflates the noise estimate) while the rank-sum test barely notices:",
            "the outlier's rank can only go up to N.")
        )),
        layout_columns(
          col_widths = c(7, 5),
          card(full_screen = TRUE, card_header("From values to ranks"),
               plotOutput(ns("rank_plot"), height = "430px"),
               axes_note("top: the actual values; bottom: their ranks (1 = smallest).", "group.",
                         "The number above each point is its rank. Ranks keep the order but forget the distances,",
                         "so an extreme value is just 'the largest'.")),
          div(uiOutput(ns("rank_boxes")))
        )
      ),

      # ==== 2. Null distribution ==============================================
      nav_panel(
        "2 · Shuffling labels", value = "null",
        explain(
          p("If the treatment does nothing, the labels 'Control' and 'Treatment' are arbitrary: any way of splitting",
            "the ranks into two groups was equally likely. So we can work out the p-value by ", strong("shuffling"),
            ": reassign the labels at random, recompute the rank sum, and repeat."),
          p("For small samples the computer can list every possible split exactly (left chart). The p-value is the",
            "share of splits whose rank sum is at least as far from the expected value as ours.",
            "Use the buttons to shuffle for yourself and watch your shuffles build up the same shape (right).")
        ),
        layout_columns(
          col_widths = c(6, 6),
          card(card_header("Every possible rank sum, if there's no difference"),
               plotOutput(ns("exact_plot"), height = "340px"),
               axes_note("rank sum of the Treatment group.", "probability of that rank sum when there's no real difference.",
                         "Black line: your rank sum. Orange bars: sums at least as extreme (their total = p-value).")),
          card(card_header("Shuffle it yourself"),
               div(class = "d-flex gap-2 mb-2",
                   actionButton(ns("shuffle1"), "Shuffle once", class = "btn-outline-primary btn-sm"),
                   actionButton(ns("shuffle100"), "Shuffle 100 times", class = "btn-outline-primary btn-sm"),
                   actionButton(ns("shuffle_reset"), "Clear", class = "btn-outline-secondary btn-sm")),
               plotOutput(ns("last_shuffle"), height = "80px"),
               plotOutput(ns("shuffle_hist"), height = "230px"),
               axes_note("rank sum after shuffling the labels.", "number of shuffles.",
                         "Top strip: the latest shuffle (ranks 1 to N, coloured by their new label).",
                         "Orange: shuffles at least as extreme as your real data."))
        )
      ),

      # ==== 3. Compare ========================================================
      nav_panel(
        "3 · t-test vs Wilcoxon", value = "compare",
        explain(
          p("Which test should you use? A good test has two properties: it rarely raises a false alarm when there",
            "is no difference (about 5% of the time), and it often detects a real difference (high ", strong("power"), ")."),
          p("This simulation runs 500 experiments with no difference and 500 with your chosen difference, and",
            "counts how often each test says 'significant'. Try each shape of noise:"),
          tags$ul(
            tags$li(strong("Bell-shaped:"), " both tests keep false alarms near 5%; the t-test has a little more power."),
            tags$li(strong("Skewed or heavy-tailed:"), " the Wilcoxon test often has noticeably more power, because",
                    "extreme values inflate the t-test's noise estimate.")
          ),
          p(class = "small text-muted",
            "Strictly, the rank-sum test asks whether one group tends to have larger values than the other,",
            "not whether the averages differ. When both groups have the same shape, that's a shift in location.")
        ),
        layout_columns(
          col_widths = c(4, 8),
          card(card_header("Shape of the noise"),
               plotOutput(ns("shape_plot"), height = "260px"),
               axes_note("value.", "how common each value is.", "A large sample from the chosen noise shape.")),
          card(card_header("How often each test says 'significant'"),
               plotOutput(ns("compare_plot"), height = "320px"),
               axes_note("scenario.", "share of the 500 experiments with p < 0.05.",
                         "Left pair: should be near 5% (dotted line). Right pair: higher is better (power)."))
        )
      ),

      # ==== 4. Signed-rank ====================================================
      nav_panel(
        "4 · Paired: signed-rank", value = "signed",
        explain(withMathJax(
          p("For paired data (the same people measured twice) the rank-based version is the ",
            strong("Wilcoxon signed-rank test"), ". Take each person's change, rank the changes by their size",
            "(ignoring the sign), then add up the ranks of the ", em("positive"), " changes."),
          math_words(
            "$$V = \\sum (\\text{ranks of } |d| \\text{ for people whose change is positive}) \\qquad \\text{expected} = \\frac{n(n+1)}{4}$$",
            "if there's no real change, each change is equally likely to be up or down, so the positive changes",
            "should collect about half of the total rank. Big changes in one direction push V away from that.",
            symbols = list("\\(d\\)" = "one person's change", "\\(|d|\\)" = "its size, ignoring the sign",
                           "\\(n\\)" = "number of people")
          )
        )),
        layout_columns(
          col_widths = c(7, 5),
          card(full_screen = TRUE, card_header("Each person's change, ordered by size"),
               plotOutput(ns("signed_plot"), height = "400px"),
               axes_note("rank of the change's size (1 = smallest change, ignoring direction).",
                         "the change (after − before).",
                         "Blue bars: increases; orange: decreases. The number on each bar is its rank.",
                         "V adds up the blue ranks.")),
          div(uiOutput(ns("signed_boxes")))
        )
      ),

      # ==== 5. Kruskal-Wallis =================================================
      nav_panel(
        "5 · More groups: Kruskal-Wallis", value = "kw",
        explain(withMathJax(
          p("The ", strong("Kruskal-Wallis test"), " is the rank-based version of one-way ANOVA. Rank all observations",
            "together, then check whether some groups have systematically higher or lower average ranks."),
          math_words(
            "$$H = \\frac{12}{N(N+1)} \\sum_j n_j \\Big(\\bar{R}_j - \\frac{N+1}{2}\\Big)^2$$",
            "for each group, measure how far its average rank is from the overall average rank (the middle,",
            "(N+1)/2), square it, and weight by group size. It's ANOVA's 'between groups' idea applied to ranks.",
            "The 12/(N(N+1)) just puts H on a standard scale (a chi-squared distribution).",
            symbols = list("\\(\\bar{R}_j\\)" = "average rank in group j", "\\(n_j\\)" = "size of group j",
                           "\\(N\\)" = "total observations")
          )
        )),
        layout_columns(
          col_widths = c(7, 5),
          card(full_screen = TRUE, card_header("Ranks by group"),
               plotOutput(ns("kw_plot"), height = "400px"),
               axes_note("group.", "rank of each observation among all of them (1 = smallest).",
                         "Black dashes: each group's average rank. Dotted line: the overall middle rank.")),
          div(uiOutput(ns("kw_boxes")))
        )
      ),

      # ==== 6. Real data ======================================================
      nav_panel(
        "6 · Real data", value = "real",
        explain(p(
          "Each dataset gets the rank-based test and its classic counterpart side by side: rank-sum vs t-test",
          "for two groups, signed-rank vs paired t-test for paired data, Kruskal-Wallis vs ANOVA for several groups.",
          "When they agree, good. When they don't, look at the plot for skew or outliers."
        )),
        layout_columns(
          col_widths = c(6, 6),
          card(card_header("The data"),
               plotOutput(ns("real_plot"), height = "380px"),
               axes_note("group.", "measured value.",
                         "Dots: observations (gray lines join paired measurements). Black: average and 95% CI.")),
          card(card_header("Results"), uiOutput(ns("real_result")))
        )
      ),

      nav_panel("Quiz", value = "quiz", quiz_panel(ns))
    )
  )
}

wilcox_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    SD <- 2

    # ---- shared two-group sample --------------------------------------------
    seed <- reactiveVal(91)
    observeEvent(input$resample, seed(sample.int(1e6, 1)))
    dat <- reactive({
      d <- sim_two_groups(input$n, input$shift, SD, seed(), input$shape)
      if (input$outlier > 0) {
        i <- which(d$group == "Treatment")
        j <- i[which.max(d$y[i])]
        d$y[j] <- d$y[j] + input$outlier
      }
      d$rank <- rank(d$y)
      d
    })
    rsum <- reactive(sum(dat()$rank[dat()$group == "Treatment"]))
    expected <- reactive(input$n * (2 * input$n + 1) / 2)
    wx <- reactive({
      d <- dat()
      suppressWarnings(wilcox.test(d$y[d$group == "Treatment"], d$y[d$group == "Control"]))
    })

    # ---- 1. ranks -----------------------------------------------------------
    output$rank_plot <- renderPlot({
      d <- dat()
      long <- rbind(data.frame(panel = "Actual values", x = d$y, group = d$group, rank = d$rank),
                    data.frame(panel = "Ranks", x = d$rank, group = d$group, rank = d$rank))
      long$panel <- factor(long$panel, levels = c("Actual values", "Ranks"))
      ggplot(long, aes(x, group, colour = group)) +
        geom_point(size = 4, alpha = 0.85) +
        ggrepel::geom_text_repel(aes(label = rank), nudge_y = 0.25, direction = "x", size = 3.3, colour = sv_col$ink,
                                 min.segment.length = Inf, seed = 1, box.padding = 0.1) +
        scale_colour_manual(values = two_cols, guide = "none") +
        facet_wrap(~panel, ncol = 1, scales = "free_x") +
        labs(x = NULL, y = NULL) +
        theme_statviz() + theme(strip.text = element_text(face = "bold", colour = sv_col$ink, hjust = 0))
    })

    output$rank_boxes <- renderUI({
      d <- dat()
      tp <- t.test(y ~ group, data = d)$p.value
      tagList(
        sv_box("Treatment rank sum", rsum(),
               p(sprintf("expected %s if there's no difference (range %d to %d)", fmt(expected(), 1),
                         input$n * (input$n + 1) / 2, input$n * (3 * input$n + 1) / 2)), tone = "blue"),
        layout_column_wrap(
          width = 1 / 2, fill = FALSE,
          sv_box("Rank-sum p-value", fmt_p(wx()$p.value), tone = if (wx()$p.value < 0.05) "green" else "gray"),
          sv_box("t-test p-value", fmt_p(tp), tone = if (tp < 0.05) "green" else "gray")
        )
      )
    })

    # ---- 2. null distribution -----------------------------------------------
    exact <- reactive({
      n <- input$n; off <- n * (n + 1) / 2
      u <- 0:(n * n)
      data.frame(R = u + off, p = dwilcox(u, n, n))
    })

    output$exact_plot <- renderPlot({
      e <- exact(); obs <- rsum(); E <- expected()
      e$extreme <- abs(e$R - E) >= abs(obs - E) - 1e-9
      pval <- sum(e$p[e$extreme])
      ggplot(e, aes(R, p, fill = extreme)) +
        geom_col(width = 0.9) +
        geom_vline(xintercept = obs, colour = sv_col$ink, linewidth = 0.9) +
        scale_fill_manual(values = c(`FALSE` = "#9ec5f4", `TRUE` = sv_col$accent), guide = "none") +
        annotate("label", x = Inf, y = Inf, hjust = 1.05, vjust = 1.3, fill = sv_col$surface, label.size = 0,
                 label = p_lab(pval), size = 4) +
        labs(x = "Treatment rank sum", y = "Probability") + theme_statviz()
    })

    shuffles <- reactiveVal(numeric(0))
    last_lab <- reactiveVal(NULL)
    observeEvent(list(input$n, input$shift, input$shape, input$outlier, input$resample, input$shuffle_reset),
                 { shuffles(numeric(0)); last_lab(NULL) })
    do_shuffle <- function(times) {
      d <- dat()
      sums <- vapply(seq_len(times), function(i) {
        lab <- sample(d$group)
        if (i == times) last_lab(lab)
        sum(d$rank[lab == "Treatment"])
      }, 0)
      shuffles(c(shuffles(), sums))
    }
    observeEvent(input$shuffle1, do_shuffle(1))
    observeEvent(input$shuffle100, do_shuffle(100))

    output$last_shuffle <- renderPlot({
      lab <- last_lab()
      validate(need(!is.null(lab), "Press a shuffle button."))
      d <- dat()
      ggplot(data.frame(rank = d$rank, g = lab), aes(rank, 1, fill = g)) +
        geom_tile(width = 0.9, height = 0.8) +
        geom_text(aes(label = rank), colour = "#ffffff", size = 3.2, fontface = "bold") +
        scale_fill_manual(values = two_cols, name = NULL) +
        labs(x = NULL, y = NULL) +
        theme_void() + theme(legend.position = "right")
    })

    output$shuffle_hist <- renderPlot({
      s <- shuffles()
      validate(need(length(s) > 0, ""))
      obs <- rsum(); E <- expected()
      df <- data.frame(R = s, extreme = abs(s - E) >= abs(obs - E) - 1e-9)
      ggplot(df, aes(R, fill = extreme)) +
        geom_histogram(binwidth = 1, colour = sv_col$surface) +
        geom_vline(xintercept = obs, colour = sv_col$ink, linewidth = 0.9) +
        scale_fill_manual(values = c(`FALSE` = "#9ec5f4", `TRUE` = sv_col$accent), guide = "none") +
        annotate("label", x = Inf, y = Inf, hjust = 1.05, vjust = 1.3, fill = sv_col$surface, label.size = 0, size = 3.6,
                 label = sprintf("%d shuffles, %.1f%% as extreme", length(s), 100 * mean(df$extreme))) +
        coord_cartesian(xlim = range(exact()$R)) +
        labs(x = "Treatment rank sum", y = "Shuffles") + theme_statviz(11)
    })

    # ---- 3. compare ---------------------------------------------------------
    output$shape_plot <- renderPlot({
      set.seed(1)
      v <- draw_noise(20000, 1, input$shape_c)
      ggplot(data.frame(v = v[abs(v) < 8]), aes(v)) +
        geom_density(fill = "#cde2fb", colour = sv_col$primary, linewidth = 1) +
        labs(x = "Value", y = "Density") + theme_statviz(11)
    })

    comp <- reactiveVal(NULL)
    observeEvent(list(input$shape_c, input$n_c, input$shift_c), comp(NULL), ignoreInit = TRUE)
    observeEvent(input$run_c, {
      n <- input$n_c; shape <- input$shape_c
      run <- function(shift, offset) {
        r <- vapply(seq_len(500), function(i) {
          d <- sim_two_groups(n, shift, 1, offset + i, shape)
          a <- d$y[d$group == "Control"]; b <- d$y[d$group == "Treatment"]
          c(t.test(b, a)$p.value < 0.05, suppressWarnings(wilcox.test(b, a)$p.value) < 0.05)
        }, numeric(2))
        rowMeans(r)
      }
      withProgress(message = "Running experiments", value = 0.2, {
        null <- run(0, 20000); incProgress(0.4)
        alt <- run(input$shift_c, 40000)
      })
      comp(data.frame(
        scenario = factor(rep(c("No real difference\n(false alarm rate)", "Real difference\n(power)"), each = 2),
                          levels = c("No real difference\n(false alarm rate)", "Real difference\n(power)")),
        test = rep(c("t-test", "Wilcoxon rank-sum"), 2),
        rate = c(null, alt)))
    })

    output$compare_plot <- renderPlot({
      df <- comp()
      validate(need(!is.null(df), "Press 'Run 500 experiments' in the sidebar."))
      ggplot(df, aes(scenario, rate, fill = test)) +
        geom_col(position = position_dodge(0.75), width = 0.7) +
        geom_text(aes(label = sprintf("%.0f%%", 100 * rate)), position = position_dodge(0.75), vjust = -0.4,
                  colour = sv_col$ink, size = 4) +
        geom_hline(yintercept = 0.05, linetype = "dotted", colour = sv_col$ink2) +
        scale_fill_manual(values = c(`t-test` = sv_col$primary, `Wilcoxon rank-sum` = sv_col$good), name = NULL) +
        scale_y_continuous(limits = c(0, 1.05), labels = function(v) paste0(100 * v, "%")) +
        labs(x = NULL, y = "Said 'significant'") + theme_statviz()
    })

    # ---- 4. signed-rank -----------------------------------------------------
    seed_p <- reactiveVal(101)
    observeEvent(input$resample_p, seed_p(sample.int(1e6, 1)))
    diffs <- reactive({ set.seed(seed_p()); input$change + draw_noise(input$np, input$sd_p, input$shape_p) })

    output$signed_plot <- renderPlot({
      d <- diffs(); r <- rank(abs(d))
      df <- data.frame(r = r, d = d, up = d > 0)
      ggplot(df, aes(r, d, fill = up)) +
        geom_hline(yintercept = 0, colour = sv_col$ink2) +
        geom_col(width = 0.8) +
        geom_text(aes(label = r, vjust = ifelse(up, -0.4, 1.3)), size = 3.4, colour = sv_col$ink) +
        scale_fill_manual(values = c(`TRUE` = sv_col$primary, `FALSE` = sv_col$accent), guide = "none") +
        scale_y_continuous(expand = expansion(mult = 0.12)) +
        scale_x_continuous(breaks = function(l) unique(round(pretty(l)))) +
        labs(x = "Rank of the size of the change", y = "Change (after − before)") + theme_statviz()
    })

    output$signed_boxes <- renderUI({
      d <- diffs(); n <- length(d)
      V <- sum(rank(abs(d))[d > 0])
      sw <- suppressWarnings(wilcox.test(d))
      pt <- t.test(d)
      tagList(
        sv_box("V: sum of ranks of increases", V,
               p(sprintf("expected %s if there's no change (range 0 to %d)", fmt(n * (n + 1) / 4, 1), n * (n + 1) / 2)),
               tone = "blue"),
        layout_column_wrap(
          width = 1 / 2, fill = FALSE,
          sv_box("Signed-rank p-value", fmt_p(sw$p.value), tone = if (sw$p.value < 0.05) "green" else "gray"),
          sv_box("Paired t-test p-value", fmt_p(pt$p.value), tone = if (pt$p.value < 0.05) "green" else "gray")
        )
      )
    })

    # ---- 5. Kruskal-Wallis --------------------------------------------------
    seed_k <- reactiveVal(111)
    observeEvent(input$resample_k, seed_k(sample.int(1e6, 1)))
    kdat <- reactive({
      d <- sim_k_groups(input$k, input$n_k, input$spread, 2, seed_k(), input$shape_k)
      d$rank <- rank(d$y)
      d
    })

    output$kw_plot <- renderPlot({
      d <- kdat(); m <- tapply(d$rank, d$group, mean); N <- nrow(d)
      set.seed(5)
      ggplot(d, aes(group, rank)) +
        geom_hline(yintercept = (N + 1) / 2, linetype = "dotted", colour = sv_col$ink2) +
        geom_jitter(width = 0.15, height = 0, colour = sv_col$primary, alpha = 0.7, size = 2.3) +
        geom_segment(data = data.frame(x = seq_along(m), m = m), aes(x = x - 0.3, xend = x + 0.3, y = m, yend = m),
                     inherit.aes = FALSE, colour = sv_col$ink, linewidth = 1.4) +
        labs(x = "Group", y = "Rank") + theme_statviz()
    })

    output$kw_boxes <- renderUI({
      d <- kdat()
      kw <- kruskal.test(y ~ group, data = d)
      av <- summary(aov(y ~ group, data = d))[[1]]$`Pr(>F)`[1]
      tagList(
        sv_box("H statistic", fmt(kw$statistic), p(sprintf("compared with a chi-squared distribution, %d df", kw$parameter)),
               tone = "blue"),
        layout_column_wrap(
          width = 1 / 2, fill = FALSE,
          sv_box("Kruskal-Wallis p", fmt_p(kw$p.value), tone = if (kw$p.value < 0.05) "green" else "gray"),
          sv_box("ANOVA p", fmt_p(av), tone = if (av < 0.05) "green" else "gray")
        ),
        p(class = "small text-muted mt-2",
          "Like ANOVA, a significant result says some group differs. Follow up with pairwise rank-sum tests",
          "using a multiple-comparison adjustment (e.g. Holm or Dunn's test).")
      )
    })

    # ---- 6. real data -------------------------------------------------------
    ds <- reactive(wx_datasets[[input$dataset]])
    output$blurb <- renderUI(p(class = "small text-muted", ds()$blurb))

    output$real_plot <- renderPlot({
      d <- ds()$data; d$g <- factor(d[[ds()$group]])
      g <- group_points_plot(d, y = ds()$y, group = "g", ylab = ds()$y, jitter = if (ds()$type == "paired") 0 else 0.12)
      if (ds()$type == "paired") g <- g + geom_line(aes(group = .data[[ds()$id]]), colour = sv_col$muted, alpha = 0.6)
      g
    })

    output$real_result <- renderUI({
      d <- ds()$data; y <- ds()$y; d$g <- factor(d[[ds()$group]]); lv <- levels(d$g)
      row <- function(name, stat, p) tags$tr(tags$td(name), tags$td(stat), tags$td(strong(fmt_p(p))))
      if (ds()$type == "k") {
        kw <- kruskal.test(d[[y]], d$g)
        av <- summary(aov(d[[y]] ~ d$g))[[1]]
        res <- list(row("Kruskal-Wallis (ranks)", sprintf("H = %s", fmt(kw$statistic)), kw$p.value),
                    row("One-way ANOVA (means)", sprintf("F = %s", fmt(av$`F value`[1])), av$`Pr(>F)`[1]))
        extra <- NULL
      } else {
        a <- d[[y]][d$g == lv[1]]; b <- d[[y]][d$g == lv[2]]
        if (ds()$type == "paired") {
          o <- order(d[[ds()$id]][d$g == lv[1]]); a <- a[o]; b <- b[order(d[[ds()$id]][d$g == lv[2]])]
          w <- suppressWarnings(wilcox.test(b, a, paired = TRUE, conf.int = TRUE, exact = FALSE))
          t <- t.test(b, a, paired = TRUE)
          res <- list(row("Wilcoxon signed-rank (ranks)", sprintf("V = %s", fmt(w$statistic, 1)), w$p.value),
                      row("Paired t-test (means)", sprintf("t = %s", fmt(t$statistic)), t$p.value))
        } else {
          w <- suppressWarnings(wilcox.test(b, a, conf.int = TRUE, exact = FALSE))
          t <- t.test(b, a)
          res <- list(row("Wilcoxon rank-sum (ranks)", sprintf("W = %s", fmt(w$statistic, 1)), w$p.value),
                      row("Welch t-test (means)", sprintf("t = %s", fmt(t$statistic)), t$p.value))
        }
        extra <- p(class = "small",
                   sprintf("Typical shift (%s − %s), estimated from ranks: %s (95%% CI %s to %s).",
                           lv[2], lv[1], fmt(w$estimate), fmt(w$conf.int[1]), fmt(w$conf.int[2])),
                   "This is the Hodges-Lehmann estimate: the median of all differences between one value from each group.")
      }
      tagList(
        tags$table(class = "table table-sm",
                   tags$thead(tags$tr(tags$th("Test"), tags$th("Statistic"), tags$th("p-value"))),
                   tags$tbody(res)),
        extra,
        p(class = "small text-muted",
          "Ties (equal values) are common in real data; the rank-based p-values then use a normal approximation.")
      )
    })

    quiz_server(input, output, session, list(
      list(q = "The Wilcoxon rank-sum test works with...",
           choices = c("the group averages", "the ranks (order) of the observations",
                       "the largest value in each group", "the standard deviations"),
           answer = 2, why = "It replaces values by their ranks, so only the ordering matters."),
      list(q = "Why is the rank-sum test less affected by an extreme outlier?",
           choices = c("It deletes outliers automatically", "An outlier's rank can be at most N, however extreme its value",
                       "It uses the median of each group", "It doesn't; it's equally affected"),
           answer = 2, why = "Ranks cap the influence of any one point: the most extreme value is simply 'the largest'."),
      list(q = "With perfectly bell-shaped data, compared with the t-test the Wilcoxon test is...",
           choices = c("much more powerful", "slightly less powerful", "invalid", "identical"),
           answer = 2, why = "It loses only a little power with normal data (about 5%), and often wins with skewed or heavy-tailed data."),
      list(q = "Which test is the rank-based version of a paired t-test?",
           choices = c("Wilcoxon rank-sum", "Wilcoxon signed-rank", "Kruskal-Wallis", "Tukey's HSD"),
           answer = 2, why = "Signed-rank ranks the sizes of within-pair changes and sums the ranks of the positive ones."),
      list(q = "Kruskal-Wallis is to one-way ANOVA as...",
           choices = c("the rank-sum test is to the t-test", "the t-test is to regression",
                       "the median is to the mode", "power is to α"),
           answer = 1, why = "Both are rank-based versions of a means-based test, for several groups and for two groups respectively.")
    ))
  })
}
