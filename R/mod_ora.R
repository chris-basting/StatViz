# Over-representation analysis (ORA) module ---------------------------------------
#
# Tabs:
#   1. Is the overlap surprising? - hypergeometric test, 2x2 table, simulation
#   2. Testing many pathways      - p-value histogram, Benjamini-Hochberg FDR
#   3. Pitfalls                   - background universe and the DE cut-off
#   4. Reading results            - enrichment dot plot from a simulated experiment
#   5. Quiz

ora_ui <- function(id) {
  ns <- NS(id)

  layout_sidebar(
    sidebar = sidebar(
      width = 300,
      show_on(ns, "idea",
        sliderInput(ns("N"), "Genes measured (the universe)", 1000, 30000, 10000, step = 500),
        sliderInput(ns("n"), "Differentially expressed (DE) genes", 10, 2000, 300, step = 10),
        sliderInput(ns("K"), "Genes in the pathway", 5, 500, 50, step = 5),
        sliderInput(ns("k"), "DE genes in the pathway (the overlap)", 0, 50, 6, step = 1)
      ),
      show_on(ns, "multiple",
        sliderInput(ns("m"), "Pathways tested", 50, 2000, 500, step = 50),
        sliderInput(ns("n_true"), "Pathways that are truly enriched", 0, 50, 10, step = 1),
        sliderInput(ns("strength"), "How enriched they are (fold)", 1, 6, 3, step = 0.25),
        sliderInput(ns("alpha"), "Significance level", 0.01, 0.2, 0.05, step = 0.01),
        actionButton(ns("resample_m"), "Run again", icon = icon("shuffle"),
                     class = "btn-outline-secondary btn-sm w-100")
      ),
      show_on(ns, "pitfalls",
        radioButtons(ns("universe"), "Background (universe)", c(
          "Genes measured in the experiment (10,000)" = "measured",
          "Every gene in the genome (20,000)" = "genome"
        )),
        sliderInput(ns("cut_pit"), "DE cut-off (adjusted p below)", 0.001, 0.2, 0.05, step = 0.001)
      ),
      show_on(ns, "results",
        sliderInput(ns("cut"), "DE cut-off (adjusted p below)", 0.001, 0.2, 0.05, step = 0.001),
        radioButtons(ns("direction"), "Which DE genes", c("Up and down" = "both", "Up only" = "up", "Down only" = "down"))
      ),
      show_on(ns, "quiz", p("Test yourself on the ideas from the other tabs."))
    ),

    navset_underline(
      id = ns("tab"),

      # ==== 1. Idea ============================================================
      nav_panel(
        "1 · Is the overlap surprising?", value = "idea",
        explain(withMathJax(
          p("A typical experiment (e.g. RNA-seq comparing treated and control cells) ends with a list of",
            strong("differentially expressed (DE) genes"), ": maybe 300 out of 10,000 measured. That list is",
            "hard to interpret gene by gene, so we ask a broader question: ",
            em("are genes from a known pathway (e.g. 'interferon response') over-represented in my list?")),
          p("If 3% of all genes are DE, we'd expect about 3% of any pathway's genes to be DE just by chance. ",
            strong("Over-representation analysis (ORA)"), " checks whether the pathway's overlap is",
            "surprisingly large. Think of the genes as marbles in a jar: ", strong("N"), " marbles, ",
            strong("K"), " of them painted (the pathway). Grab ", strong("n"), " at random (the DE list). How",
            "likely is it to get ", strong("k"), " or more painted ones?"),
          math_words(
            "$$P(X \\ge k) = \\sum_{i \\ge k} \\frac{\\binom{K}{i} \\binom{N-K}{n-i}}{\\binom{N}{n}}$$",
            "count the ways to grab i painted and n − i unpainted marbles, divide by all the ways to grab n, and",
            "add up over every i at least as big as what we saw. This is the ", strong("hypergeometric test"),
            ", the same as a one-sided ", strong("Fisher's exact test"), " on the 2×2 table below.",
            symbols = list("\\(\\binom{a}{b}\\)" = "'a choose b': ways to pick b items out of a",
                           "\\(N\\)" = "genes measured", "\\(K\\)" = "genes in the pathway",
                           "\\(n\\)" = "DE genes", "\\(k\\)" = "DE genes in the pathway")
          ),
          math_words(
            "$$\\text{fold enrichment} = \\frac{k}{\\text{expected}} = \\frac{k}{K \\cdot n / N}$$",
            "how many times more DE genes the pathway has than chance would give. A fold of 1 means 'what",
            "we'd expect'; 4 means four times as many."
          )
        )),
        layout_columns(
          col_widths = c(4, 8),
          div(
            card(card_header("The 2 × 2 table"), plotOutput(ns("tab2x2"), height = "230px"),
                 card_footer(class = "small sv-axes", "Rows: in the pathway or not. Columns: DE or not.",
                             "Every gene measured lands in exactly one cell.")),
            uiOutput(ns("ora_boxes"))
          ),
          card(full_screen = TRUE, card_header("How big would the overlap be by chance?"),
               plotOutput(ns("hyper_plot"), height = "380px"),
               axes_note("number of the pathway's genes that land in the DE list.",
                         "probability, if the DE list were a random pick of genes.",
                         "Black line: your overlap. Orange bars: overlaps at least as large (their total is the p-value).",
                         "Dotted line: the expected overlap."))
        )
      ),

      # ==== 2. Multiple testing ================================================
      nav_panel(
        "2 · Testing many pathways", value = "multiple",
        explain(withMathJax(
          p("Pathway databases contain hundreds or thousands of gene sets, and we test them all. At",
            "p < 0.05, about 5% of the pathways with ", em("no"), " real enrichment will look significant anyway.",
            "With 500 pathways, that's about 25 false discoveries."),
          p("The standard fix is to control the ", strong("false discovery rate (FDR)"),
            ": among the pathways we call significant, keep the expected share of false ones below e.g. 5%.",
            "The ", strong("Benjamini-Hochberg"), " procedure does this:"),
          math_words(
            "$$\\text{sort the p-values: } p_{(1)} \\le p_{(2)} \\le \\dots \\le p_{(m)}; \\quad \\text{call significant all } i \\le i^* \\text{, the largest } i \\text{ with } p_{(i)} \\le \\frac{i}{m}\\,\\alpha$$",
            "line the p-values up from smallest to largest and compare each with a threshold that rises along",
            "the line (the diagonal in the right-hand chart). Everything up to the last p-value under the line",
            "is called significant. The adjusted p-value reported by software (the ", strong("q-value"),
            ") is the smallest FDR at which each pathway would be called significant.",
            symbols = list("\\(m\\)" = "number of pathways tested", "\\(\\alpha\\)" = "target FDR (e.g. 0.05)",
                           "\\(p_{(i)}\\)" = "the i-th smallest p-value")
          )
        )),
        layout_columns(
          col_widths = c(6, 6),
          card(card_header("p-values of all pathways"),
               plotOutput(ns("p_hist"), height = "330px"),
               axes_note("p-value.", "number of pathways.",
                         "Pathways with no real enrichment spread evenly from 0 to 1 (flat part). True enrichment",
                         "piles up near 0. Colour: what's really going on (we know because we simulated it).")),
          card(card_header("Benjamini-Hochberg, step by step"),
               plotOutput(ns("bh_plot"), height = "330px"),
               axes_note("rank of the p-value (1 = smallest), first 60 shown.", "p-value.",
                         "Diagonal: the rising BH threshold i/m × α. Dashed: the fixed α cut-off.",
                         "Everything left of the vertical line is called significant by BH."))
        ),
        uiOutput(ns("fdr_table"))
      ),

      # ==== 3. Pitfalls ========================================================
      nav_panel(
        "3 · Pitfalls", value = "pitfalls",
        explain(
          p(strong("1. The background (universe) matters."), " The question is 'compared with what?'. The right",
            "universe is the genes that ", em("could"), " have been called DE: those actually measured and expressed in",
            "your experiment. Using every gene in the genome inflates enrichment for pathways whose genes are",
            "mostly expressed in your tissue (e.g. metabolic pathways in liver), because the many silent genes in",
            "the genome make the expected overlap look smaller than it should be. Switch the universe in the",
            "sidebar: those pathways jump above the diagonal, while others can drop below it. Either way, the",
            "p-values no longer mean what they should."),
          p(strong("2. The DE cut-off is arbitrary, and ORA only sees genes that pass it."), " A pathway where",
            em("many"), " genes shift a little (e.g. every cell-cycle gene up by 20%) may have few genes passing any",
            "cut-off, so ORA misses it. The right-hand chart follows the five pathways that really changed in our",
            "simulated experiment. Three of them are broad and subtle and never become significant; that's",
            "the gap GSEA (next page) fills.")
        ),
        layout_columns(
          col_widths = c(6, 6),
          card(full_screen = TRUE, card_header("Same data, two backgrounds"),
               plotOutput(ns("universe_plot"), height = "380px"),
               axes_note("−log10 p with the correct background (genes measured). Higher = more significant.",
                         "−log10 p with the background chosen in the sidebar.",
                         "Each dot is a pathway. On the diagonal = no change. Above = looks more significant than it should.",
                         "Labelled: pathways whose genes are mostly expressed in this tissue.")),
          card(full_screen = TRUE, card_header("What different cut-offs find"),
               plotOutput(ns("cutoff_plot"), height = "380px"),
               axes_note("DE cut-off (adjusted p, log scale).", "−log10 q for the pathway (higher = more significant).",
                         "Each line is one of the five pathways that truly changed. Dashed line: q = 0.05.",
                         "Orange line: your cut-off."))
        )
      ),

      # ==== 4. Results =========================================================
      nav_panel(
        "4 · Reading results", value = "results",
        explain(
          p("Here is a full ORA on a simulated experiment: 10,000 genes, 40 pathways. ORA results are usually",
            "shown as a ", strong("dot plot"), " like the one below. For each pathway:"),
          tags$ul(
            tags$li(strong("Fold enrichment"), " (or 'gene ratio'): how over-represented the pathway is."),
            tags$li(strong("Count:"), " how many DE genes it contains (dot size)."),
            tags$li(strong("Significance:"), " FDR-adjusted p-value (colour).")
          ),
          p("Look for large, dark dots far to the right. A huge fold enrichment from 2 genes in a tiny pathway is",
            "much weaker evidence than a moderate fold from 30 genes. Also try 'Up only' and 'Down only': mixing",
            "directions can hide pathways that go clearly one way.")
        ),
        layout_columns(
          col_widths = c(7, 5),
          card(full_screen = TRUE, card_header("Top pathways"),
               plotOutput(ns("dotplot"), height = "470px"),
               axes_note("fold enrichment (DE genes found ÷ expected by chance).", "pathway (top 15 by p-value).",
                         "Dot size: number of DE genes in the pathway. Colour: −log10 q (darker = more significant).",
                         "Dashed line: fold 1 (no enrichment).")),
          card(card_header("Results table"), uiOutput(ns("res_table")))
        )
      ),

      nav_panel("Quiz", value = "quiz", quiz_panel(ns))
    )
  )
}

ora_server <- function(id) {
  moduleServer(id, function(input, output, session) {

    # ---- 1. idea ------------------------------------------------------------
    observe({
      mx <- min(input$K, input$n)
      updateSliderInput(session, "k", max = mx, value = min(isolate(input$k), mx))
    })
    vals <- reactive({
      N <- input$N; n <- min(input$n, N); K <- min(input$K, N); k <- min(input$k, K, n)
      list(N = N, n = n, K = K, k = k, exp = K * n / N, p = ora_test(k, K, n, N))
    })

    output$tab2x2 <- renderPlot({
      v <- vals()
      df <- data.frame(row = factor(c("In pathway", "In pathway", "Not in pathway", "Not in pathway"),
                                    levels = c("Not in pathway", "In pathway")),
                       col = factor(c("DE", "Not DE", "DE", "Not DE"), levels = c("DE", "Not DE")),
                       n = c(v$k, v$K - v$k, v$n - v$k, v$N - v$K - v$n + v$k),
                       key = c(TRUE, FALSE, FALSE, FALSE))
      ggplot(df, aes(col, row)) +
        geom_tile(aes(fill = key), colour = sv_col$surface, linewidth = 2) +
        geom_text(aes(label = format(n, big.mark = ",")), size = 6, fontface = "bold", colour = sv_col$ink) +
        scale_fill_manual(values = c(`TRUE` = "#cde2fb", `FALSE` = sv_col$neutral), guide = "none") +
        scale_x_discrete(position = "top") + labs(x = NULL, y = NULL) +
        theme_statviz() + theme(panel.grid = element_blank(), axis.text = element_text(size = 12))
    })

    output$ora_boxes <- renderUI({
      v <- vals()
      layout_column_wrap(
        width = 1 / 2, fill = FALSE,
        sv_box("Expected overlap", fmt(v$exp, 1), p(sprintf("observed %d", v$k)), tone = "gray"),
        sv_box("Fold enrichment", if (v$exp > 0) fmt(v$k / v$exp, 1) else "–", tone = "blue"),
        sv_box("p-value", fmt_p(v$p), p(if (v$p < 0.05) "a surprising overlap" else "could easily be chance"),
               tone = if (v$p < 0.05) "green" else "orange")
      )
    })

    output$hyper_plot <- renderPlot({
      v <- vals()
      hi <- max(v$k + 3, qhyper(0.9999, v$K, v$N - v$K, v$n))
      x <- 0:min(hi, v$K, v$n)
      df <- data.frame(x = x, p = dhyper(x, v$K, v$N - v$K, v$n), tail = x >= v$k)
      ggplot(df, aes(x, p, fill = tail)) +
        geom_col(width = 0.8) +
        geom_vline(xintercept = v$k, colour = sv_col$ink, linewidth = 0.9) +
        geom_vline(xintercept = v$exp, colour = sv_col$ink2, linetype = "dotted") +
        scale_fill_manual(values = c(`FALSE` = "#9ec5f4", `TRUE` = sv_col$accent), guide = "none") +
        annotate("label", x = Inf, y = Inf, hjust = 1.05, vjust = 1.3, fill = sv_col$surface, label.size = 0, size = 4,
                 label = sprintf("P(overlap ≥ %d) = %s", v$k, fmt_p(v$p))) +
        labs(x = "Overlap (DE genes in the pathway)", y = "Probability") + theme_statviz()
    })

    # ---- 2. multiple testing ------------------------------------------------
    seed_m <- reactiveVal(161)
    observeEvent(input$resample_m, seed_m(sample.int(1e6, 1)))
    multi <- reactive({
      set.seed(seed_m())
      m <- input$m; nt <- min(input$n_true, m); N <- 10000; n <- 300
      K <- round(exp(runif(m, log(15), log(300))))
      truth <- seq_len(m) <= nt
      k <- ifelse(truth, rbinom(m, K, pmin(0.9, input$strength * n / N)), rhyper(m, K, N - K, n))
      p <- ora_test(k, K, n, N)
      data.frame(K = K, k = k, p = p, q = p.adjust(p, "BH"), truth = truth)
    })
    truth_cols <- c(`Truly enriched` = sv_col$primary, `No real enrichment` = "#b8b7b2")

    output$p_hist <- renderPlot({
      d <- multi()
      d$what <- factor(ifelse(d$truth, "Truly enriched", "No real enrichment"), levels = names(truth_cols))
      ggplot(d, aes(p, fill = what)) +
        geom_histogram(breaks = seq(0, 1, 0.05), colour = sv_col$surface) +
        geom_vline(xintercept = input$alpha, linetype = "22", colour = sv_col$ink) +
        scale_fill_manual(values = truth_cols, name = NULL) +
        labs(x = "p-value", y = "Pathways") + theme_statviz()
    })

    output$bh_plot <- renderPlot({
      d <- multi(); m <- nrow(d); a <- input$alpha
      d <- d[order(d$p), ]; d$rank <- seq_len(m)
      below <- which(d$p <= d$rank / m * a); istar <- if (length(below)) max(below) else 0
      show <- head(d, 60)
      show$what <- factor(ifelse(show$truth, "Truly enriched", "No real enrichment"), levels = names(truth_cols))
      ggplot(show, aes(rank, p)) +
        geom_abline(intercept = 0, slope = a / m, colour = sv_col$good, linewidth = 1) +
        geom_hline(yintercept = a, linetype = "22", colour = sv_col$ink2) +
        { if (istar > 0) geom_vline(xintercept = istar + 0.5, colour = sv_col$accent, linewidth = 0.9) } +
        geom_point(aes(colour = what), size = 2.5) +
        scale_colour_manual(values = truth_cols, name = NULL) +
        coord_cartesian(ylim = c(0, min(1, max(a * 1.6, a * 60 / m * 1.3)))) +
        labs(x = "Rank", y = "p-value") + theme_statviz()
    })

    output$fdr_table <- renderUI({
      d <- multi(); a <- input$alpha
      row <- function(name, sig) {
        tp <- sum(sig & d$truth); fp <- sum(sig & !d$truth)
        tags$tr(tags$td(name), tags$td(sum(sig)), tags$td(tp), tags$td(strong(fp)),
                tags$td(if (sum(sig)) sprintf("%.0f%%", 100 * fp / sum(sig)) else "–"),
                tags$td(sprintf("%d of %d", tp, sum(d$truth))))
      }
      card(
        card_header("What each approach calls significant"),
        tags$table(class = "table table-sm",
                   tags$thead(tags$tr(tags$th("Rule"), tags$th("Called significant"), tags$th("Truly enriched"),
                                      tags$th("False discoveries"), tags$th("Share false"), tags$th("Found"))),
                   tags$tbody(row(sprintf("Raw p < %s", a), d$p < a),
                              row(sprintf("BH q < %s (FDR control)", a), d$q < a),
                              row(sprintf("Bonferroni p < %s / %d", a, nrow(d)), d$p < a / nrow(d)))),
        p(class = "small text-muted mb-0",
          "Bonferroni controls the chance of even one false discovery, which is very strict, so it usually finds fewer",
          "real pathways. BH accepts a small, controlled share of false ones in exchange for more discoveries.")
      )
    })

    # ---- shared simulated experiment (tabs 3-4) -----------------------------
    exper <- sim_experiment(seed = 1)
    sets <- attr(exper, "sets"); active <- attr(exper, "active")

    # genome-wide version of the sets: add members that aren't expressed in this tissue
    genome_sets <- local({
      set.seed(11)
      tissue <- c("Oxidative phosphorylation", "Fatty acid metabolism", "Peroxisome", "Bile acid metabolism",
                  "Xenobiotic metabolism", "Cholesterol homeostasis", "Heme metabolism", "Adipogenesis")
      unexpressed <- sprintf("SILENT%05d", 1:10000)
      lapply(setNames(names(sets), names(sets)), function(s) {
        rate <- if (s %in% tissue) 0.95 else runif(1, 0.3, 0.6)     # share of members expressed here
        extra <- round(length(sets[[s]]) * (1 / rate - 1))
        c(sets[[s]], sample(unexpressed, extra))
      })
    })
    tissue_sets <- c("Oxidative phosphorylation", "Fatty acid metabolism", "Peroxisome", "Bile acid metabolism",
                     "Xenobiotic metabolism", "Cholesterol homeostasis", "Heme metabolism", "Adipogenesis")

    de_at <- function(cut, direction = "both") {
      sel <- exper$padj < cut
      if (direction == "up") sel <- sel & exper$t > 0
      if (direction == "down") sel <- sel & exper$t < 0
      exper$gene[sel]
    }

    output$universe_plot <- renderPlot({
      de <- de_at(input$cut_pit)
      validate(need(length(de) > 0, "No genes pass this cut-off."))
      a <- ora_all(de, exper$gene, sets)
      b <- if (input$universe == "genome") ora_all(de, c(exper$gene, sprintf("SILENT%05d", 1:10000)), genome_sets) else a
      df <- merge(a[, c("set", "p")], b[, c("set", "p")], by = "set", suffixes = c("_ok", "_sel"))
      df$x <- -log10(df$p_ok); df$y <- -log10(df$p_sel)
      df$tissue <- df$set %in% tissue_sets
      lim <- c(0, max(c(df$x, df$y, 3)) * 1.05)
      ggplot(df, aes(x, y)) +
        geom_abline(colour = sv_col$ink2, linetype = "22") +
        geom_point(aes(colour = tissue), size = 3, alpha = 0.85) +
        ggrepel::geom_text_repel(data = df[df$tissue & df$y > df$x + 0.3, ], aes(label = set), size = 3.3,
                                 colour = sv_col$ink, seed = 1, max.overlaps = 20) +
        scale_colour_manual(values = c(`TRUE` = sv_col$accent, `FALSE` = sv_col$primary),
                            labels = c(`TRUE` = "mostly expressed in this tissue", `FALSE` = "other pathways"), name = NULL) +
        coord_equal(xlim = lim, ylim = lim) +
        labs(x = "−log10 p, correct background", y = "−log10 p, chosen background") + theme_statviz()
    })

    output$cutoff_plot <- renderPlot({
      cuts <- 10^seq(log10(0.001), log10(0.2), length.out = 12)
      df <- do.call(rbind, lapply(cuts, function(c) {
        de <- de_at(c); r <- ora_all(de, exper$gene, sets)
        data.frame(cut = c, set = active$set, q = r$q[match(active$set, r$set)])
      }))
      df$lq <- -log10(df$q)
      lab <- df[df$cut == max(df$cut), ]
      cols <- setNames(sv_groups[1:5], active$set)
      ggplot(df, aes(cut, lq, colour = set)) +
        geom_hline(yintercept = -log10(0.05), linetype = "22", colour = sv_col$ink2) +
        geom_vline(xintercept = input$cut_pit, colour = sv_col$accent, linewidth = 0.9) +
        geom_line(linewidth = 1.1) + geom_point(size = 1.8) +
        ggrepel::geom_text_repel(data = lab, aes(label = set), hjust = 0, nudge_x = 0.15, direction = "y",
                                 size = 3.2, seed = 1, segment.colour = NA) +
        scale_x_log10(breaks = c(0.001, 0.01, 0.1), labels = c("0.001", "0.01", "0.1"),
                      expand = expansion(mult = c(0.02, 0.9))) +
        scale_colour_manual(values = cols, guide = "none") +
        labs(x = "DE cut-off (adjusted p)", y = "−log10 q") + theme_statviz()
    })

    # ---- 4. results ---------------------------------------------------------
    res <- reactive({
      de <- de_at(input$cut, input$direction)
      validate(need(length(de) > 0, "No genes pass this cut-off."))
      ora_all(de, exper$gene, sets)
    })

    output$dotplot <- renderPlot({
      r <- head(res(), 15)
      r$set <- factor(r$set, levels = rev(r$set))
      r$lq <- -log10(r$q)
      ggplot(r, aes(fold, set)) +
        geom_vline(xintercept = 1, linetype = "22", colour = sv_col$ink2) +
        geom_point(aes(size = k, colour = lq)) +
        scale_colour_gradientn(colours = seq_blues[2:7], name = "−log10 q") +
        scale_size_area(max_size = 10, name = "DE genes") +
        labs(x = "Fold enrichment", y = NULL) + theme_statviz() + theme(legend.position = "right")
    })

    output$res_table <- renderUI({
      r <- head(res(), 12)
      n_de <- length(de_at(input$cut, input$direction))
      tagList(
        p(class = "small", strong(n_de), " DE genes at this cut-off (out of 10,000)."),
        tags$table(class = "table table-sm small",
                   tags$thead(tags$tr(tags$th("Pathway"), tags$th("k / K"), tags$th("Fold"), tags$th("q"))),
                   tags$tbody(lapply(seq_len(nrow(r)), function(i) tags$tr(
                     tags$td(if (r$q[i] < 0.05) strong(r$set[i]) else r$set[i]),
                     tags$td(sprintf("%d / %d", r$k[i], r$K[i])), tags$td(fmt(r$fold[i], 1)),
                     tags$td(fmt_p(r$q[i])))))),
        p(class = "small text-muted mb-0",
          "Bold: q < 0.05. Truly changed in this simulation: ",
          paste(active$set, collapse = ", "), ".")
      )
    })

    quiz_server(input, output, session, list(
      list(q = "ORA asks whether...",
           choices = c("a pathway's genes are over-represented among the DE genes, compared with chance",
                       "every gene in the pathway changed", "the pathway is biologically important",
                       "the DE genes are correlated"),
           answer = 1, why = "It compares the observed overlap with what a random gene list would give."),
      list(q = "10,000 genes measured, 500 DE, a pathway of 40 genes. The expected overlap by chance is...",
           choices = c("0.5", "2", "20", "40"),
           answer = 2, why = "Expected = K × n / N = 40 × 500 / 10,000 = 2."),
      list(q = "You test 1,000 pathways with no real enrichment at p < 0.05. Roughly how many look significant?",
           choices = c("0", "5", "50", "500"),
           answer = 3, why = "5% of 1,000 = 50 false positives, which is why FDR correction (e.g. Benjamini-Hochberg) is needed."),
      list(q = "Which background (universe) should ORA use?",
           choices = c("Every gene in the genome", "Only the DE genes", "The genes measured/expressed in the experiment",
                       "Only genes in the pathway"),
           answer = 3, why = "Only genes that could have been called DE belong in the comparison; otherwise tissue-specific pathways look falsely enriched."),
      list(q = "A pathway where most genes shift slightly in the same direction is often missed by ORA because...",
           choices = c("ORA only looks at genes that pass the DE cut-off", "ORA ignores pathway size",
                       "ORA requires at least 1,000 genes", "small shifts are never biologically relevant"),
           answer = 1, why = "Few genes pass the cut-off individually, so the coordinated shift is invisible. GSEA uses the whole ranked list instead.")
    ))
  })
}
