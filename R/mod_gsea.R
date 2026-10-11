# Gene set enrichment analysis (GSEA) module ----------------------------------------
#
# Tabs:
#   1. Rank every gene        - no cut-off; where do a pathway's genes fall?
#   2. The running sum        - enrichment score step by step
#   3. Is it significant?     - permutation null, NES, p-value
#   4. ORA vs GSEA            - broad-subtle vs sparse-strong pathways
#   5. Reading results        - NES chart and the classic enrichment plot
#   6. Quiz

gsea_ui <- function(id) {
  ns <- NS(id)

  layout_sidebar(
    sidebar = sidebar(
      width = 300,
      show_on(ns, c("rank", "es", "perm"),
        h6("Simulated experiment"),
        sliderInput(ns("N"), "Genes measured", 500, 10000, 2000, step = 500),
        sliderInput(ns("K"), "Genes in the pathway", 10, 300, 50, step = 5),
        sliderInput(ns("shift"), "Shift of the pathway's genes (in t units)", -2, 2, 0.6, step = 0.1),
        sliderInput(ns("frac"), "Share of the pathway's genes that shift", 0.1, 1, 0.8, step = 0.05),
        actionButton(ns("resample"), "Run the experiment again", icon = icon("shuffle"),
                     class = "btn-outline-secondary btn-sm w-100")
      ),
      show_on(ns, "es",
        hr(),
        sliderInput(ns("pos"), "Walk down the list to position", 1, 10000, 10000, step = 1),
        radioButtons(ns("weight"), "Step size for each pathway gene", c(
          "Equal steps (classic Kolmogorov-Smirnov, p = 0)" = 0,
          "Proportional to the gene's statistic (standard GSEA, p = 1)" = 1
        ), selected = 1)
      ),
      show_on(ns, "perm",
        hr(),
        sliderInput(ns("n_perm"), "Number of permutations", 100, 5000, 1000, step = 100),
        actionButton(ns("run_perm"), "Run permutations", icon = icon("play"), class = "btn-primary w-100")
      ),
      show_on(ns, "compare",
        sliderInput(ns("broad"), "Broad pathway: shift of each of its 80 genes", 0, 1.5, 0.5, step = 0.05),
        sliderInput(ns("sparse"), "Sparse pathway: shift of 10 of its 80 genes", 0, 8, 5, step = 0.25),
        sliderInput(ns("cut_c"), "ORA cut-off (adjusted p below)", 0.001, 0.2, 0.05, step = 0.001),
        actionButton(ns("resample_c"), "Run the experiment again", icon = icon("shuffle"),
                     class = "btn-outline-secondary btn-sm w-100")
      ),
      show_on(ns, "results",
        selectInput(ns("set_sel"), "Show the enrichment plot for", pathway_names, selected = "Oxidative phosphorylation")
      ),
      show_on(ns, "quiz", p("Test yourself on the ideas from the other tabs."))
    ),

    navset_underline(
      id = ns("tab"),

      # ==== 1. Ranking =========================================================
      nav_panel(
        "1 · Rank every gene", value = "rank",
        explain(
          p("ORA needs a cut-off to decide which genes are 'DE', and everything below it is thrown away.",
            strong("Gene set enrichment analysis (GSEA)"), " uses ", em("every"), " gene instead:"),
          tags$ol(
            tags$li("Rank all genes by how much they changed, from most increased at the top to most decreased at the",
                    "bottom. The ranking statistic is usually the t-statistic or a signed fold change."),
            tags$li("Ask whether a pathway's genes are scattered randomly through the list, or bunched toward the top",
                    "or the bottom.")
          ),
          p("Even if no single gene in the pathway is dramatically changed, many genes shifting a little in the same",
            "direction will crowd toward one end of the list. Try a small shift (e.g. 0.3) with most of the pathway",
            "affected, and the bunching is still visible.")
        ),
        card(full_screen = TRUE, card_header("All genes, ranked"),
             plotOutput(ns("rank_metric"), height = "260px"),
             plotOutput(ns("rank_barcode"), height = "90px"),
             axes_note("rank (1 = most increased, right end = most decreased).",
                       "top: the ranking statistic (t) for each gene; bottom: a tick for each gene in the pathway.",
                       "If the pathway didn't matter, its ticks would be spread evenly along the whole list.")),
        uiOutput(ns("rank_boxes"))
      ),

      # ==== 2. Running sum =====================================================
      nav_panel(
        "2 · The running sum", value = "es",
        explain(withMathJax(
          p("To turn 'bunched toward one end' into a number, GSEA walks down the ranked list keeping a running",
            "score. Each pathway gene ", strong("steps it up"), "; every other gene ", strong("steps it down"),
            " a little. The steps are sized so that the walk always ends back at 0."),
          math_words(
            "$$\\text{RS}(i) = \\sum_{\\substack{\\text{pathway genes}\\\\ \\text{up to } i}} \\frac{|t_j|^p}{\\sum_{\\text{pathway}} |t|^p} \\; - \\sum_{\\substack{\\text{other genes}\\\\ \\text{up to } i}} \\frac{1}{N - K}$$",
            "pathway genes add their share of the pathway's total 'weight'; other genes each subtract the same",
            "small amount. If the pathway's genes are concentrated near the top, the score climbs early and",
            "high. The ", strong("enrichment score (ES)"), " is the furthest the walk gets from 0, up or down.",
            "With p = 1 (standard), strongly changed genes take bigger steps; with p = 0 every gene counts equally.",
            symbols = list("\\(t_j\\)" = "ranking statistic of gene j", "\\(N\\)" = "genes in the list",
                           "\\(K\\)" = "genes in the pathway", "\\(p\\)" = "weight exponent")
          ),
          p("The pathway genes at or before the peak form the ", strong("leading edge"),
            ": the subset driving the enrichment, often the most interesting genes to follow up.",
            "Use the position slider to walk down the list yourself.")
        )),
        card(full_screen = TRUE, card_header("The GSEA enrichment plot"),
             plotOutput(ns("es_curve"), height = "280px"),
             plotOutput(ns("es_barcode"), height = "70px"),
             plotOutput(ns("es_metric"), height = "150px"),
             axes_note("rank in the list (1 = most increased).",
                       "top: running enrichment score; middle: pathway genes; bottom: ranking statistic.",
                       "Green dashed line and dot: the enrichment score (the peak). Shaded ticks: the leading edge.",
                       "Orange line: your position in the walk.")),
        uiOutput(ns("es_boxes"))
      ),

      # ==== 3. Permutation =====================================================
      nav_panel(
        "3 · Is it significant?", value = "perm",
        explain(withMathJax(
          p("Is our ES large? We compare it with ES values from random 'pathways' of the same size: pick K genes",
            "at random, compute their ES, and repeat many times. That's the ", strong("null distribution"),
            ". The p-value is the share of random pathways with an ES at least as extreme (in the same direction)."),
          math_words(
            "$$\\text{NES} = \\frac{\\text{ES}}{\\text{average } |\\text{ES}| \\text{ of random sets (same sign)}}$$",
            "the ", strong("normalized enrichment score"), " rescales ES so that pathways of different sizes can be",
            "compared: big pathways have smaller random ES values, so their raw ES isn't comparable. NES around ±1 is",
            "typical of chance; |NES| above about 1.5 is usually interesting."
          ),
          p(class = "small text-muted",
            "Here we shuffle genes ('gene-set permutation', as in preranked GSEA). With enough samples (roughly 7+",
            "per group), the original GSEA shuffles the ", em("sample labels"), " instead, which also respects the",
            "fact that genes in a pathway tend to be correlated.")
        )),
        layout_columns(
          col_widths = c(7, 5),
          card(card_header("Your ES compared with random pathways"),
               plotOutput(ns("perm_hist"), height = "360px"),
               axes_note("enrichment score of a random gene set of the same size.", "number of random sets.",
                         "Black line: your pathway's ES. Orange: random sets at least as extreme in the same direction.",
                         "Random ES values are never near 0: ES is the biggest swing of the walk, up or down, so random sets",
                         "land on either side, and each side is compared separately.")),
          div(uiOutput(ns("perm_boxes")))
        )
      ),

      # ==== 4. Compare =========================================================
      nav_panel(
        "4 · ORA vs GSEA", value = "compare",
        explain(
          p("The two methods see different things. Here is one experiment with two pathways of 80 genes each:"),
          tags$ul(
            tags$li(strong("Broad:"), " every gene shifts a little. Few genes pass a DE cut-off, so ORA struggles,",
                    "but GSEA sees all 80 genes crowding toward the top."),
            tags$li(strong("Sparse:"), " 10 genes change a lot; the rest do nothing. Those 10 pass the cut-off, so ORA",
                    "picks it up. GSEA usually does too, though a handful of genes in a big set gives a weaker signal.")
          ),
          p("Rules of thumb: use GSEA when you have a ranking for all genes (e.g. from DESeq2 or limma). Use ORA when",
            "all you have is a gene list (e.g. genes near GWAS hits, or a cluster of co-expressed genes).")
        ),
        layout_columns(
          col_widths = c(6, 6),
          card(card_header("Broad pathway (small shifts, all genes)"),
               plotOutput(ns("cmp_broad"), height = "220px"), uiOutput(ns("cmp_broad_res"))),
          card(card_header("Sparse pathway (big shifts, few genes)"),
               plotOutput(ns("cmp_sparse"), height = "220px"), uiOutput(ns("cmp_sparse_res")))
        ),
        p(class = "small sv-axes mx-1",
          strong("Charts:"), " running enrichment score along the ranked list of 5,000 genes, with the pathway's genes as",
          "ticks below. Orange ticks: genes that pass the ORA cut-off.")
      ),

      # ==== 5. Results =========================================================
      nav_panel(
        "5 · Reading results", value = "results",
        explain(
          p("Full GSEA on the same simulated experiment as the ORA page (10,000 genes, 40 pathways, 1,000",
            "permutations each). Pathways with positive NES are enriched among increased genes; negative NES among",
            "decreased genes. GSEA results are often filtered at FDR q < 0.25, a deliberately loose threshold",
            "recommended by its authors for exploration."),
          p("Five pathways were really changed in this simulation: interferon alpha response and EMT (few genes,",
            "strong), and E2F targets, oxidative phosphorylation and hypoxia (many genes, subtle). Compare with the",
            "ORA page, which only found the first two.")
        ),
        layout_columns(
          col_widths = c(6, 6),
          card(full_screen = TRUE, card_header("Normalized enrichment score for every pathway"),
               plotOutput(ns("nes_plot"), height = "560px"),
               axes_note("NES: positive = enriched among increased genes, negative = among decreased genes.", "pathway.",
                         "Strong colour: FDR q < 0.05; pale: not significant. Red: up; blue: down.")),
          div(
            card(full_screen = TRUE, card_header(textOutput(ns("sel_title"), inline = TRUE)),
                 plotOutput(ns("sel_curve"), height = "200px"),
                 plotOutput(ns("sel_barcode"), height = "60px"),
                 axes_note("rank (1 = most increased).", "running enrichment score; ticks = pathway genes.")),
            card(card_header("Top results"), uiOutput(ns("gsea_table")))
          )
        )
      ),

      nav_panel("Quiz", value = "quiz", quiz_panel(ns))
    )
  )
}

gsea_server <- function(id) {
  moduleServer(id, function(input, output, session) {

    # ---- plot builders ------------------------------------------------------
    curve_plot <- function(rs, es, peak, cursor = NULL, compact = FALSE) {
      df <- data.frame(i = seq_along(rs), rs = rs)
      g <- ggplot(df, aes(i, rs)) +
        geom_hline(yintercept = 0, colour = sv_col$ink2) +
        geom_hline(yintercept = es, colour = sv_col$good, linetype = "22") +
        geom_line(colour = sv_col$primary, linewidth = 1) +
        annotate("point", x = peak, y = es, colour = sv_col$good, size = 3.5)
      if (!is.null(cursor)) g <- g + geom_vline(xintercept = cursor, colour = sv_col$accent, linewidth = 0.9) +
        annotate("point", x = cursor, y = rs[cursor], colour = sv_col$accent, size = 3.5)
      g + scale_x_continuous(expand = c(0, 0)) +
        labs(x = NULL, y = if (compact) "ES" else "Running enrichment score") +
        theme_statviz(if (compact) 10 else 12) + theme(axis.text.x = element_blank())
    }
    barcode_plot <- function(N, pos, lead = NULL, highlight = NULL) {
      df <- data.frame(x = pos, lead = if (is.null(lead)) FALSE else pos %in% lead,
                       hl = if (is.null(highlight)) FALSE else pos %in% highlight)
      col <- ifelse(df$hl, sv_col$accent, ifelse(df$lead, sv_col$ink, "#8a8984"))
      ggplot(df) +
        geom_segment(aes(x = x, xend = x, y = 0, yend = 1), colour = col, linewidth = 0.4) +
        scale_x_continuous(limits = c(1, N), expand = c(0, 0)) +
        theme_void() + theme(plot.background = element_rect(fill = sv_col$surface, colour = NA))
    }
    metric_plot <- function(metric) {
      df <- data.frame(i = seq_along(metric), m = metric)
      ggplot(df, aes(i, m)) +
        geom_area(aes(fill = m > 0), alpha = 0.7) +
        geom_hline(yintercept = 0, colour = sv_col$ink2) +
        scale_fill_manual(values = c(`TRUE` = sv_col$pos, `FALSE` = sv_col$neg), guide = "none") +
        scale_x_continuous(expand = c(0, 0)) +
        labs(x = "Rank in list", y = "Statistic (t)") + theme_statviz(11)
    }

    # ---- simulated ranked list (tabs 1-3) -----------------------------------
    seed <- reactiveVal(171)
    observeEvent(input$resample, seed(sample.int(1e6, 1)))
    sim <- reactive({
      set.seed(seed())
      N <- input$N; K <- min(input$K, N - 1)
      members <- sample(N, K)
      hit <- sample(members, round(input$frac * K))
      t <- rnorm(N); t[hit] <- t[hit] + input$shift
      o <- order(t, decreasing = TRUE)
      list(metric = t[o], pos = sort(match(members, o)), N = N, K = K)
    })
    observeEvent(input$N, updateSliderInput(session, "pos", max = input$N, value = input$N))
    w_exp <- reactive(as.numeric(input$weight %||% 1))
    walk <- reactive({
      s <- sim(); hits <- seq_len(s$N) %in% s$pos
      running_sum(s$metric, hits, w_exp())
    })
    lead_pos <- function(pos, es, peak) if (es >= 0) pos[pos <= peak] else pos[pos >= peak]

    output$rank_metric <- renderPlot(metric_plot(sim()$metric))
    output$rank_barcode <- renderPlot(barcode_plot(sim()$N, sim()$pos))
    output$rank_boxes <- renderUI({
      s <- sim()
      top <- mean(s$pos <= s$N / 2)
      layout_column_wrap(
        width = 1 / 3, fill = FALSE,
        sv_box("Pathway genes in the top half", sprintf("%.0f%%", 100 * top), p("50% expected by chance"),
               tone = if (abs(top - 0.5) > 0.15) "green" else "gray"),
        sv_box("Average rank of pathway genes", format(round(mean(s$pos)), big.mark = ","),
               p(sprintf("%s expected by chance", format(round((s$N + 1) / 2), big.mark = ","))), tone = "blue"),
        sv_box("Pathway genes with |t| > 3", sum(abs(s$metric[s$pos]) > 3),
               p("roughly what passes a strict DE cut-off"), tone = "gray")
      )
    })

    # ---- 2. running sum -----------------------------------------------------
    output$es_curve <- renderPlot({
      w <- walk(); cur <- min(input$pos, length(w$rs))
      curve_plot(w$rs, w$es, w$peak, cursor = cur)
    })
    output$es_barcode <- renderPlot({
      s <- sim(); w <- walk()
      barcode_plot(s$N, s$pos, lead = lead_pos(s$pos, w$es, w$peak))
    })
    output$es_metric <- renderPlot(metric_plot(sim()$metric) + geom_vline(xintercept = min(input$pos, sim()$N),
                                                                           colour = sv_col$accent, linewidth = 0.9))
    output$es_boxes <- renderUI({
      s <- sim(); w <- walk(); cur <- min(input$pos, s$N)
      layout_column_wrap(
        width = 1 / 3, fill = FALSE,
        sv_box("Enrichment score (ES)", fmt(w$es), p(sprintf("peak at rank %s", format(w$peak, big.mark = ","))),
               tone = "green"),
        sv_box("Leading edge", length(lead_pos(s$pos, w$es, w$peak)), p(sprintf("of the %d pathway genes", s$K)),
               tone = "blue"),
        sv_box("Running score here", fmt(w$rs[cur]), p(sprintf("after %s genes, %d of them in the pathway",
                                                               format(cur, big.mark = ","), sum(s$pos <= cur))),
               tone = "orange")
      )
    })

    # ---- 3. permutations ----------------------------------------------------
    perm <- reactiveVal(NULL)
    observeEvent(list(input$N, input$K, input$shift, input$frac, input$resample, input$weight), perm(NULL), ignoreInit = TRUE)
    observeEvent(input$run_perm, {
      s <- sim(); hits <- seq_len(s$N) %in% s$pos
      withProgress(message = "Shuffling genes", value = 0.3, {
        perm(gsea_one(s$metric, hits, n_perm = input$n_perm, weight = w_exp(), seed = seed()))
      })
    })

    output$perm_hist <- renderPlot({
      r <- perm()
      validate(need(!is.null(r), "Press 'Run permutations' in the sidebar."))
      df <- data.frame(es = r$null)
      df$extreme <- if (r$es >= 0) df$es >= r$es else df$es <= r$es
      ggplot(df, aes(es, fill = extreme)) +
        geom_histogram(bins = 50, colour = sv_col$surface) +
        geom_vline(xintercept = r$es, colour = sv_col$ink, linewidth = 1) +
        scale_fill_manual(values = c(`FALSE` = "#9ec5f4", `TRUE` = sv_col$accent), guide = "none") +
        labs(x = "Enrichment score of random sets", y = "Random sets") + theme_statviz()
    })

    output$perm_boxes <- renderUI({
      r <- perm()
      if (is.null(r)) return(p(class = "text-muted", "Results appear here after running the permutations."))
      tagList(
        sv_box("Enrichment score", fmt(r$es), tone = "blue"),
        sv_box("Normalized (NES)", fmt(r$nes), p("|NES| > 1.5 is usually interesting"), tone = "blue"),
        sv_box("p-value", fmt_p(r$p),
               p(sprintf("%d permutations; the smallest possible p is about %s", length(r$null), fmt(1 / (length(r$null) / 2 + 1), 4))),
               tone = if (r$p < 0.05) "green" else "orange")
      )
    })

    # ---- 4. compare ---------------------------------------------------------
    seed_c <- reactiveVal(181)
    observeEvent(input$resample_c, seed_c(sample.int(1e6, 1)))
    cmp <- reactive({
      set.seed(seed_c())
      N <- 5000; genes <- seq_len(N)
      broad <- sample(genes, 80); sparse <- sample(setdiff(genes, broad), 80)
      t <- rnorm(N)
      t[broad] <- t[broad] + input$broad
      strong <- sample(sparse, 10); t[strong] <- t[strong] + input$sparse
      # a background of strongly changed genes, as in a real experiment
      bg <- sample(setdiff(genes, c(broad, sparse)), 100); t[bg] <- t[bg] + sample(c(-4.5, 4.5), 100, TRUE)
      padj <- p.adjust(2 * pnorm(-abs(t)), "BH")
      de <- genes[padj < input$cut_c]
      o <- order(t, decreasing = TRUE); metric <- t[o]
      one <- function(m) {
        pos <- sort(match(m, o)); hits <- genes %in% pos
        g <- gsea_one(metric, hits, n_perm = 1000, seed = 1)
        rs <- running_sum(metric, hits, 1)
        list(pos = pos, rs = rs, gsea = g, de_pos = sort(match(intersect(m, de), o)),
             ora_k = length(intersect(m, de)), ora_p = ora_test(length(intersect(m, de)), 80, length(de), N))
      }
      list(N = N, n_de = length(de), broad = one(broad), sparse = one(sparse))
    })

    cmp_plot <- function(x, N) {
      top <- curve_plot(x$rs$rs, x$rs$es, x$rs$peak, compact = TRUE)
      ticks <- data.frame(x = x$pos, de = x$pos %in% x$de_pos)
      ylim <- range(c(x$rs$rs, 0))
      h <- diff(ylim) * 0.18
      top + geom_segment(data = ticks, aes(x = x, xend = x, y = ylim[1] - 1.6 * h, yend = ylim[1] - 0.4 * h,
                                           colour = de), inherit.aes = FALSE, linewidth = 0.5) +
        scale_colour_manual(values = c(`TRUE` = sv_col$accent, `FALSE` = "#8a8984"), guide = "none")
    }
    cmp_res <- function(x, n_de) {
      tags$table(class = "table table-sm small mb-0",
                 tags$tbody(
                   tags$tr(tags$td("ORA"), tags$td(sprintf("%d of 80 genes pass the cut-off (%d DE overall)", x$ora_k, n_de)),
                           tags$td(strong(p_lab(x$ora_p)))),
                   tags$tr(tags$td("GSEA"), tags$td(sprintf("NES %s", fmt(x$gsea$nes))), tags$td(strong(p_lab(x$gsea$p))))))
    }
    output$cmp_broad <- renderPlot(cmp_plot(cmp()$broad, cmp()$N))
    output$cmp_sparse <- renderPlot(cmp_plot(cmp()$sparse, cmp()$N))
    output$cmp_broad_res <- renderUI(cmp_res(cmp()$broad, cmp()$n_de))
    output$cmp_sparse_res <- renderUI(cmp_res(cmp()$sparse, cmp()$n_de))

    # ---- 5. results on the shared experiment --------------------------------
    exper <- sim_experiment(seed = 1)
    sets <- attr(exper, "sets")
    full <- reactive(gsea_all(exper[, c("gene", "t")], sets, n_perm = 1000, seed = 2)) |> bindCache("gsea-full")
    ranked <- local({ o <- order(exper$t, decreasing = TRUE); list(metric = exper$t[o], genes = exper$gene[o]) })

    output$nes_plot <- renderPlot({
      r <- full()
      r$set <- factor(r$set, levels = r$set[order(r$nes)])
      r$dir <- ifelse(r$nes > 0, "up", "down")
      r$sig <- r$q < 0.05
      ggplot(r, aes(nes, set, fill = interaction(dir, sig))) +
        geom_col(width = 0.75) +
        geom_vline(xintercept = 0, colour = sv_col$ink2) +
        scale_fill_manual(values = c(up.TRUE = sv_col$pos, down.TRUE = sv_col$neg, up.FALSE = "#f3c4c3", down.FALSE = "#cde2fb"),
                          guide = "none") +
        labs(x = "NES", y = NULL) + theme_statviz(11)
    })

    sel <- reactive({
      pos <- sort(match(sets[[input$set_sel]], ranked$genes))
      hits <- seq_along(ranked$metric) %in% pos
      list(pos = pos, w = running_sum(ranked$metric, hits, 1))
    })
    output$sel_title <- renderText(sprintf("Enrichment plot: %s", input$set_sel))
    output$sel_curve <- renderPlot({ s <- sel(); curve_plot(s$w$rs, s$w$es, s$w$peak, compact = TRUE) })
    output$sel_barcode <- renderPlot({
      s <- sel(); barcode_plot(length(ranked$metric), s$pos, lead = lead_pos(s$pos, s$w$es, s$w$peak))
    })

    output$gsea_table <- renderUI({
      r <- head(full(), 10)
      tagList(
        tags$table(class = "table table-sm small",
                   tags$thead(tags$tr(tags$th("Pathway"), tags$th("Size"), tags$th("NES"), tags$th("q"),
                                      tags$th("Leading edge"))),
                   tags$tbody(lapply(seq_len(nrow(r)), function(i) tags$tr(
                     tags$td(if (r$q[i] < 0.05) strong(r$set[i]) else r$set[i]), tags$td(r$size[i]),
                     tags$td(fmt(r$nes[i])), tags$td(fmt_p(r$q[i])),
                     tags$td(sprintf("%d genes", r$leading_edge[i])))))),
        p(class = "small text-muted mb-0",
          "q-values here are Benjamini-Hochberg adjusted permutation p-values (the GSEA software uses a",
          "closely related NES-based FDR).")
      )
    })

    quiz_server(input, output, session, list(
      list(q = "What does GSEA use that ORA doesn't?",
           choices = c("A cut-off for DE genes", "The ranking of all genes, without a cut-off",
                       "Only the top 10 genes", "Gene names"),
           answer = 2, why = "GSEA ranks every gene and asks where the pathway's genes fall in that ranking."),
      list(q = "In the running sum, a pathway gene makes the score...",
           choices = c("go up", "go down", "reset to zero", "stay the same"),
           answer = 1, why = "Pathway genes step the score up; all other genes step it down slightly."),
      list(q = "A pathway has NES = −2.3 and q = 0.001. This means...",
           choices = c("the pathway is enriched among increased genes", "the pathway is enriched among decreased genes",
                       "the pathway has 2.3 times more genes than expected", "the result is not significant"),
           answer = 2, why = "Negative NES = genes concentrated at the bottom (decreased) end of the ranked list."),
      list(q = "Why is ES normalized into NES?",
           choices = c("To make it positive", "So pathways of different sizes can be compared",
                       "To remove the leading edge", "To convert it to a p-value"),
           answer = 2, why = "Random ES values depend on set size; dividing by the average random |ES| puts all sets on one scale."),
      list(q = "Which pathway is ORA most likely to miss but GSEA to find?",
           choices = c("6 genes changed hugely, the rest unchanged",
                       "Most of its 80 genes shifted slightly in the same direction",
                       "No genes changed", "All genes changed in random directions"),
           answer = 2, why = "Small coordinated shifts rarely pass a DE cut-off, but they push the whole set toward one end of the ranking.")
    ))
  })
}
