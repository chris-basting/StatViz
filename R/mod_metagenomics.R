# Metagenomics module ----------------------------------------------------------
#
# Tabs:
#   1. The big picture   - pipeline; shotgun sampling and coverage
#   2. Reads & quality   - Phred scores, trimming
#   3. Alignment         - dynamic programming (global/local), read mapping
#   4. Assembly          - k-mers, de Bruijn graphs, contigs, repeats, errors
#   5. Binning           - GC + coverage clustering into MAGs
#   6. Annotation        - gene finding and homology search
#   7. Who's there?      - taxonomic profiling, abundance bias, diversity
#   8. Quiz

mg_genomes <- data.frame(name = c("Species A", "Species B", "Species C"), kb = c(20, 10, 5))

metagenomics_ui <- function(id) {
  ns <- NS(id)

  layout_sidebar(
    sidebar = sidebar(
      width = 300,
      show_on(ns, "overview",
        h6("Shotgun sequencing a community"),
        p(class = "small text-muted", "Three species with genomes of 20, 10 and 5 kb (real ones are ~2-10 Mb)."),
        sliderInput(ns("abA"), "Species A: share of cells (%)", 1, 100, 60),
        sliderInput(ns("abB"), "Species B: share of cells (%)", 1, 100, 30),
        sliderInput(ns("abC"), "Species C: share of cells (%)", 1, 100, 10),
        sliderInput(ns("n_reads"), "Number of reads", 50, 3000, 1000, step = 50),
        sliderInput(ns("read_len"), "Read length (bases)", 50, 300, 150, step = 10),
        actionButton(ns("resample_o"), "Sequence again", icon = icon("shuffle"),
                     class = "btn-outline-secondary btn-sm w-100")
      ),
      show_on(ns, "reads",
        sliderInput(ns("q_thresh"), "Trim when quality drops below (Q)", 2, 38, 20),
        sliderInput(ns("min_len"), "Discard reads shorter than", 20, 140, 50, step = 5),
        helpText("Quality is checked in a sliding window of 4 bases, from the start of the read.")
      ),
      show_on(ns, "align",
        textInput(ns("seq_a"), "Sequence 1", "GATTACAGATTC"),
        textInput(ns("seq_b"), "Sequence 2", "GCATTACGATC"),
        radioButtons(ns("align_type"), "Alignment type", c(
          "Global (whole sequences, Needleman-Wunsch)" = "global",
          "Local (best matching part, Smith-Waterman)" = "local"
        )),
        sliderInput(ns("s_match"), "Score for a match", 1, 5, 2),
        sliderInput(ns("s_mis"), "Score for a mismatch", -5, 0, -1),
        sliderInput(ns("s_gap"), "Score for a gap", -6, 0, -2),
        helpText("Letters other than A, C, G, T are fine too (e.g. protein sequences). Up to 25 characters each.")
      ),
      show_on(ns, "assembly",
        h6("The genome"),
        sliderInput(ns("g_len"), "Genome length", 40, 300, 120, step = 10),
        sliderInput(ns("rep_len"), "Length of a repeated segment (0 = none)", 0, 30, 0),
        actionButton(ns("new_genome"), "New random genome", icon = icon("dna"),
                     class = "btn-outline-secondary btn-sm w-100 mb-3"),
        h6("Sequencing"),
        sliderInput(ns("a_read_len"), "Read length", 10, 60, 30),
        sliderInput(ns("a_cov"), "Coverage (average reads per base)", 1, 40, 12),
        sliderInput(ns("a_err"), "Sequencing error rate (%)", 0, 5, 0, step = 0.5),
        h6("Assembler"),
        sliderInput(ns("k"), "k (k-mer length)", 3, 31, 11, step = 2),
        sliderInput(ns("min_count"), "Ignore k-mers seen fewer than", 1, 5, 1)
      ),
      show_on(ns, "binning",
        sliderInput(ns("n_gen"), "Genomes in the sample", 2, 4, 3),
        sliderInput(ns("similar"), "How similar the genomes are (GC and abundance)", 0, 1, 0, step = 0.05),
        sliderInput(ns("contig_len"), "Typical contig length (bases)", 500, 10000, 3000, step = 500),
        actionButton(ns("resample_b"), "New assembly", icon = icon("shuffle"),
                     class = "btn-outline-secondary btn-sm w-100")
      ),
      show_on(ns, "annotation",
        sliderInput(ns("diverge"), "How different our organism's genes are from the database (%)", 0, 50, 10),
        sliderInput(ns("min_orf"), "Shortest gene to report (codons)", 20, 120, 50, step = 5),
        actionButton(ns("resample_n"), "New contig", icon = icon("shuffle"),
                     class = "btn-outline-secondary btn-sm w-100 mb-3"),
        uiOutput(ns("orf_pick"))
      ),
      show_on(ns, "taxonomy",
        h6("True community (share of cells, %)"),
        uiOutput(ns("tax_sliders")),
        sliderInput(ns("t_reads"), "Reads sequenced", 1000, 100000, 20000, step = 1000),
        hr(),
        radioButtons(ns("kmer_case"), "Example read", c(
          "From a species in the database" = "clear",
          "From a region shared with a relative" = "shared",
          "From an organism not in the database" = "novel"
        ))
      ),
      show_on(ns, "quiz", p("Test yourself on the ideas from the other tabs."))
    ),

    navset_underline(
      id = ns("tab"),

      # ==== 1. Overview ========================================================
      nav_panel(
        "1 · The big picture", value = "overview",
        explain(withMathJax(
          p(strong("Metagenomics"), " studies all the DNA in a sample (soil, seawater, a gut) at once,",
            "without growing each microbe in the lab first. Most microbes can't be cultured, so this is often",
            "the only way to see who's there and what they can do."),
          p("In ", strong("shotgun sequencing"), ", the DNA from every cell is broken into millions of short random",
            "fragments, and each fragment is read by a sequencer. The result is a huge pile of short ",
            strong("reads"), " (typically 150 letters of A, C, G, T) from every genome, all mixed together, with no",
            "label saying which organism each came from. Everything else is about making sense of that pile:"),
          tags$ul(
            tags$li(strong("Assembly-based route:"), " stitch overlapping reads back into longer sequences (contigs), group",
                    "contigs into genomes (binning) and find the genes on them (annotation)."),
            tags$li(strong("Read-based route:"), " compare each read directly to reference databases (alignment or k-mer",
                    "matching) to estimate who is there and in what proportions.")
          ),
          p("Abundant species and species with big genomes contribute more reads. The ", strong("coverage"),
            " of a genome is how many reads, on average, overlap each of its bases. Low coverage leaves gaps."),
          math_words(
            "$$c = \\frac{N \\times L}{G} \\qquad \\text{fraction of genome seen} \\approx 1 - e^{-c}$$",
            "coverage is (number of reads from that genome × read length) ÷ genome length. Because reads land at",
            "random, some bases are missed by chance; at coverage 1 about 37% of the genome is never seen, at",
            "coverage 5 less than 1%. (This is the Lander-Waterman model.)",
            symbols = list("\\(N\\)" = "reads from the genome", "\\(L\\)" = "read length", "\\(G\\)" = "genome length",
                           "\\(e^{-c}\\)" = "chance a given base gets no reads")
          )
        )),
        card(card_header("The metagenomics pipeline"),
             plotOutput(ns("pipeline"), height = "250px"),
             card_footer(class = "small sv-axes", "Numbers refer to the tabs on this page.")),
        layout_columns(
          col_widths = c(8, 4),
          card(full_screen = TRUE, card_header("Where the reads landed on each genome"),
               plotOutput(ns("coverage_plot"), height = "380px"),
               axes_note("position along the genome (kb).", "coverage: number of reads overlapping that position.",
                         "Orange ticks along the bottom: stretches no read covered (gaps in what we know).",
                         "Dashed line: average coverage.")),
          card(card_header("Coverage per species"), uiOutput(ns("coverage_table")))
        )
      ),

      # ==== 2. Reads & quality =================================================
      nav_panel(
        "2 · Reads & quality", value = "reads",
        explain(withMathJax(
          p("Sequencers aren't perfect. For every base they report a confidence score, the ",
            strong("Phred quality score"), " Q."),
          math_words(
            "$$Q = -10 \\log_{10}(p_{\\text{error}}) \\qquad p_{\\text{error}} = 10^{-Q/10}$$",
            "every 10 points of Q means ten times fewer errors: Q10 = 1 error in 10 bases, Q20 = 1 in 100,",
            "Q30 = 1 in 1,000. Q30 is a common bar for 'good'.",
            symbols = list("\\(p_{\\text{error}}\\)" = "chance the base call is wrong")
          ),
          p("Quality typically drops toward the end of each read (the chemistry wears out over cycles), and a few",
            "reads are poor throughout. ", strong("Quality control (QC)"), " trims the bad ends and throws away",
            "reads that end up too short. Other QC steps not shown here: removing leftover adapter sequence,",
            "and removing host DNA (e.g. human reads in a gut sample).")
        )),
        layout_columns(
          col_widths = c(7, 5),
          card(full_screen = TRUE, card_header("Quality along the read (300 simulated reads)"),
               plotOutput(ns("qual_plot"), height = "380px"),
               axes_note("position in the read (base 1 = first sequenced).", "Phred quality score Q.",
                         "Line: median quality at each position. Dark band: middle 50% of reads; light band: middle 90%.",
                         "Orange line: your trimming threshold.")),
          div(
            uiOutput(ns("qc_boxes")),
            card(card_header("Read lengths after trimming"),
                 plotOutput(ns("trim_hist"), height = "220px"),
                 axes_note("read length after trimming.", "number of reads.",
                           "Orange: reads discarded for being shorter than your minimum."))
          )
        )
      ),

      # ==== 3. Alignment =======================================================
      nav_panel(
        "3 · Alignment", value = "align",
        explain(withMathJax(
          p(strong("Alignment"), " lines up two sequences to show where they match. It's how we compare a read",
            "with a reference genome, a gene with a database, or two related genes. Sequences change by",
            strong("substitutions"), " (one letter swapped) and ", strong("insertions/deletions"),
            " (letters added or lost), so an alignment can include gaps (–)."),
          p("To find the best alignment, give each column a score (match, mismatch or gap) and find the",
            "arrangement with the highest total. Trying every arrangement is impossible, but ",
            strong("dynamic programming"), " solves it by filling in a table one cell at a time:"),
          math_words(
            "$$S_{i,j} = \\max \\begin{cases} S_{i-1,j-1} + \\text{match or mismatch} & \\text{(line up the two letters)} \\\\ S_{i-1,j} + \\text{gap} & \\text{(gap in sequence 2)} \\\\ S_{i,j-1} + \\text{gap} & \\text{(gap in sequence 1)} \\end{cases}$$",
            "each cell holds the best score for aligning the first i letters of one sequence with the first j",
            "letters of the other. It only needs the three neighbouring cells already filled in, so the whole",
            "table is built row by row. Following the arrows back from the end (the ", strong("traceback"),
            ") reads off the alignment. A ", strong("local"), " alignment also allows 0 ('start fresh here'),",
            "which finds the best-matching stretch instead of forcing the whole sequences to align.",
            symbols = list("\\(S_{i,j}\\)" = "best score so far for the first i and j letters")
          ),
          p(strong("Mapping millions of reads: "), "filling a table the size of (read × genome) for every read",
            "would take far too long. Read mappers (BWA, Bowtie2, minimap2) first ", strong("index"),
            " the reference: a lookup table of where every short word (a 'seed', e.g. 20 letters) occurs.",
            "Each read's seeds point straight to a few candidate locations, and the careful dynamic-programming",
            "alignment is only done there. This 'seed and extend' trick is also how BLAST and DIAMOND",
            "search huge databases.")
        )),
        layout_columns(
          col_widths = c(7, 5),
          card(full_screen = TRUE, card_header("The dynamic-programming table"),
               plotOutput(ns("dp_plot"), height = "470px"),
               axes_note("letters of sequence 2.", "letters of sequence 1.",
                         "Each cell: best score for aligning the sequences up to that point (darker = higher).",
                         "Outlined cells: the traceback path. Diagonal steps line up two letters; horizontal or",
                         "vertical steps are gaps.")),
          div(
            card(card_header("The alignment"), uiOutput(ns("aln_text"))),
            uiOutput(ns("aln_boxes"))
          )
        )
      ),

      # ==== 4. Assembly ========================================================
      nav_panel(
        "4 · Assembly", value = "assembly",
        explain(withMathJax(
          p(strong("Assembly"), " rebuilds long sequences from short reads, like reconstructing a shredded",
            "document from overlapping strips. Comparing every read with every other is too slow for millions of",
            "reads, so modern assemblers (MEGAHIT, metaSPAdes) use a ", strong("de Bruijn graph"), ":"),
          tags$ol(
            tags$li("Chop every read into overlapping words of length k, called ", strong("k-mers"),
                    ". The read ACGTTG with k = 4 gives ACGT, CGTT, GTTG."),
            tags$li("Make each k-mer a link between its first k−1 letters and its last k−1 letters (ACG → CGT).",
                    "k-mers that overlap by k−1 letters automatically chain together."),
            tags$li("Walk along chains with no forks. Each unbranched chain becomes a ", strong("contig"),
                    " (contiguous sequence).")
          ),
          p("Things that break contigs, which you can try in the sidebar:"),
          tags$ul(
            tags$li(strong("Repeats"), " longer than k: the graph can't tell which copy leads where, so it forks.",
                    "Larger k resolves longer repeats, but needs more coverage, because every k-mer must be seen."),
            tags$li(strong("Low coverage:"), " missing k-mers leave gaps."),
            tags$li(strong("Sequencing errors:"), " each error creates up to k false k-mers that are usually seen",
                    "only once. Ignoring rare k-mers cleans most of them up.")
          ),
          math_words(
            "$$\\text{N50} = \\text{the contig length } L \\text{ such that contigs} \\ge L \\text{ hold half of all assembled bases}$$",
            "a standard measure of how 'unbroken' an assembly is. Bigger is better."
          ),
          p(class = "small text-muted",
            "Simplifications: reads come from one strand only, and there's one genome. Real metagenome assemblers",
            "also handle both DNA strands, uneven coverage across species, and closely related strains.")
        )),
        layout_columns(
          col_widths = c(6, 6),
          card(full_screen = TRUE, card_header("The reads, piled on the genome"),
               plotOutput(ns("pileup"), height = "260px"),
               axes_note("position in the genome.", "reads, stacked so they don't overlap.",
                         "Each bar is one read where it really came from (the assembler doesn't know this!).",
                         "Orange marks: sequencing errors.")),
          card(full_screen = TRUE, card_header("How often each k-mer was seen"),
               plotOutput(ns("kmer_hist"), height = "260px"),
               axes_note("times a k-mer appears across all reads.", "number of distinct k-mers.",
                         "Blue: k-mers that really are in the genome. Orange: k-mers created by sequencing errors.",
                         "Dashed line: your 'ignore rarer than' cut-off."))
        ),
        layout_columns(
          col_widths = c(7, 5),
          card(full_screen = TRUE, card_header("The assembly graph"),
               plotOutput(ns("asm_graph"), height = "420px"),
               axes_note("(layout only: position has no meaning)", "",
                         "Each circle is a contig (an unbranched chain of k-mers), sized and labelled by its length.",
                         "Arrows: contigs that could follow each other. Forks are where the assembler had to stop.")),
          div(
            uiOutput(ns("asm_boxes")),
            card(card_header("Where the contigs belong"),
                 plotOutput(ns("asm_map"), height = "200px"),
                 axes_note("position in the true genome.", "contigs (longest at the top).",
                           "Green: the true genome. Blue: contigs placed where they match exactly.",
                           "Contigs containing errors can't be placed and aren't shown."))
          )
        )
      ),

      # ==== 5. Binning =========================================================
      nav_panel(
        "5 · Binning", value = "binning",
        explain(
          p("A metagenome assembly is a jumble of contigs from many organisms. ", strong("Binning"),
            " sorts them into groups that each (ideally) represent one genome, called a ",
            strong("metagenome-assembled genome (MAG)"), "."),
          p("Binners (MetaBAT, CONCOCT, MaxBin) use two kinds of clues that are shared by all pieces of one genome:"),
          tags$ul(
            tags$li(strong("Composition:"), " each genome has a characteristic DNA 'accent', e.g. its ",
                    strong("GC content"), " (share of G and C letters) and how often each 4-letter word appears."),
            tags$li(strong("Coverage:"), " all contigs of one genome come from the same number of cells, so they're",
                    "covered by roughly the same number of reads. Using several samples makes this even stronger:",
                    "contigs of one genome rise and fall together across samples.")
          ),
          p("Here, a clustering algorithm (k-means) groups contigs by GC content and coverage. Make the genomes",
            "more similar, or the contigs shorter (short contigs give noisier GC estimates), and bins start to",
            "mix. Bin quality is judged by ", strong("completeness"), " (how much of the genome made it in) and ",
            strong("contamination"), " (how much belongs to other genomes). In real data both are estimated",
            "from single-copy marker genes (CheckM); here we know the truth.")
        ),
        layout_columns(
          col_widths = c(6, 6),
          card(full_screen = TRUE, card_header("Contigs coloured by their true genome"),
               plotOutput(ns("bin_truth"), height = "380px"),
               axes_note("GC content of the contig.", "coverage (reads per base, log scale).",
                         "Each dot is a contig; size = length. Colour and shape = which genome it really came from.")),
          card(full_screen = TRUE, card_header("Contigs coloured by the bin they were put in"),
               plotOutput(ns("bin_pred"), height = "380px"),
               axes_note("GC content of the contig.", "coverage (log scale).",
                         "Colour and shape = bin chosen by k-means. Compare with the left: mismatches are binning errors."))
        ),
        card(card_header("How good are the bins?"), uiOutput(ns("bin_table")))
      ),

      # ==== 6. Annotation ======================================================
      nav_panel(
        "6 · Annotation", value = "annotation",
        explain(withMathJax(
          p(strong("Annotation"), " finds the genes on a contig and guesses what they do. Two steps:"),
          p(strong("1. Gene prediction."), " A protein-coding gene is an ", strong("open reading frame (ORF)"),
            ": it starts with the start codon ATG and runs in steps of three letters (", strong("codons"),
            ") until a stop codon (TAA, TAG or TGA). Each codon codes for one amino acid. DNA can be read in",
            "three frames on each strand, so the gene finder checks all of them (we show the three forward",
            "frames). Short ORFs pop up by chance in any DNA, so a minimum length is required. Real tools",
            "(Prodigal) also look for ribosome binding sites and codon usage patterns."),
          p(strong("2. Functional annotation."), " Each predicted protein is compared with databases of",
            "proteins whose function is known (UniProt, KEGG, COG, Pfam) using alignment tools such as BLAST,",
            "DIAMOND or HMMER. A strong match to a known protein suggests a similar function (",
            strong("homology"), "). Proteins with no good match are labelled 'hypothetical protein', which is",
            "often a large share of a metagenome."),
          math_words(
            "$$\\text{identity} = \\frac{\\text{matching amino acids}}{\\text{aligned length}}$$",
            "two random proteins share only about 5-10% identity. Above roughly 30-40% over a long stretch, the",
            "proteins are very likely related. Raise the divergence slider: matches weaken, and mutations can",
            "even create early stop codons that break a gene in two."
          )
        )),
        card(full_screen = TRUE, card_header("Genes found on a 2,000-base contig"),
             plotOutput(ns("orf_map"), height = "250px"),
             axes_note("position on the contig.", "reading frame (1, 2 or 3).",
                       "Arrows: open reading frames longer than your minimum. Blue = matched a known protein;",
                       "gray = no match ('hypothetical protein'). Green bands below: where the real genes are.")),
        layout_columns(
          col_widths = c(7, 5),
          card(card_header("Best database match for each gene"), uiOutput(ns("annot_table"))),
          card(card_header("From DNA to protein"), uiOutput(ns("codon_view")))
        )
      ),

      # ==== 7. Taxonomy ========================================================
      nav_panel(
        "7 · Who's there?", value = "taxonomy",
        explain(withMathJax(
          p(strong("Taxonomic profiling"), " estimates which organisms are present and in what proportions.",
            "Fast classifiers (Kraken2) break each read into k-mers and look each one up in a database built",
            "from reference genomes. Each k-mer points to the taxon (or group of taxa) whose genomes contain it."),
          tags$ul(
            tags$li("If the k-mers point to one species, the read is assigned to that species."),
            tags$li("If they're shared by several related species, the read goes to their ",
                    strong("lowest common ancestor"), " (e.g. the genus), the most specific label that's still safe."),
            tags$li("If few k-mers match anything, the read is ", strong("unclassified"),
                    ". Organisms missing from the database are invisible this way.")
          ),
          p(strong("Watch out:"), " the share of ", em("reads"), " isn't the share of ", em("cells"),
            ". A species with a genome twice as long sheds twice as many reads per cell. Dividing read counts by",
            "genome size corrects for this."),
          math_words(
            "$$\\hat{p}_i = \\frac{r_i / G_i}{\\sum_j r_j / G_j} \\qquad H = -\\sum_i p_i \\ln p_i$$",
            "estimated share of cells = reads ÷ genome size, rescaled to add to 1. The ",
            strong("Shannon diversity"), " H summarizes the community: 0 when one species dominates, higher",
            "when many species are evenly mixed.",
            symbols = list("\\(r_i\\)" = "reads assigned to taxon i", "\\(G_i\\)" = "its genome size",
                           "\\(p_i\\)" = "its share of the community")
          )
        )),
        card(card_header("Classifying one read by its k-mers"),
             plotOutput(ns("kmer_read"), height = "150px"),
             uiOutput(ns("kmer_verdict")),
             axes_note("the read, cut into overlapping k-mers (one tile each).", "",
                       "Tile colour = what the database says each k-mer belongs to.")),
        layout_columns(
          col_widths = c(8, 4),
          card(full_screen = TRUE, card_header("Truth vs estimates"),
               plotOutput(ns("tax_bars"), height = "380px"),
               axes_note("taxon.", "share of the community (%).",
                         "Green: true share of cells. Orange: share of reads (biased toward big genomes).",
                         "Blue: estimate after correcting for genome size. The unknown species can't be estimated: its reads go",
                         "unclassified or to its closest relative in the database (P. copri here). And because it's invisible,",
                         "every other species' estimated share is inflated.")),
          div(uiOutput(ns("tax_boxes")))
        )
      ),

      nav_panel("Quiz", value = "quiz", quiz_panel(ns))
    )
  )
}

metagenomics_server <- function(id) {
  moduleServer(id, function(input, output, session) {

    # ---- 1. overview --------------------------------------------------------
    output$pipeline <- renderPlot({
      box <- data.frame(
        x = c(1, 2, 3, 4, 5.5, 6.7, 7.9, 5.5, 7.3),
        y = c(0, 0, 0, 0, 0.6, 0.6, 0.6, -0.6, -0.6),
        lab = c("Sample\n(soil, gut...)", "Extract\nall DNA", "Shotgun\nsequencing", "Reads +\nQC  (2)",
                "Assembly\n(4)", "Binning\ninto MAGs (5)", "Annotation\n(6)", "Alignment /\nk-mers (3, 7)",
                "Who's there &\nwhat they do (7)"),
        route = c(rep("shared", 4), rep("assembly", 3), rep("read", 2)))
      arr <- data.frame(x = c(1, 2, 3, 4, 4, 5.5, 6.7, 5.5), y = c(0, 0, 0, 0, 0, 0.6, 0.6, -0.6),
                        xend = c(2, 3, 4, 5.5, 5.5, 6.7, 7.9, 7.3), yend = c(0, 0, 0, 0.6, -0.6, 0.6, 0.6, -0.6))
      arr <- transform(arr, x = x + 0.42 * sign(xend - x), xend = xend - 0.42 * sign(xend - x))
      ggplot() +
        geom_segment(data = arr, aes(x, y, xend = xend, yend = yend), colour = sv_col$ink2,
                     arrow = arrow(length = unit(0.18, "cm"), type = "closed"), linewidth = 0.6) +
        geom_label(data = box, aes(x, y, label = lab, fill = route), colour = sv_col$ink, size = 3.6,
                   label.padding = unit(0.45, "lines"), label.size = 0, lineheight = 0.95) +
        scale_fill_manual(values = c(shared = "#f0efec", assembly = "#cde2fb", read = "#e6f6ef"), guide = "none") +
        annotate("text", x = 6.7, y = 1.05, label = "assembly-based route", colour = sv_col$primary, fontface = "bold", size = 3.6) +
        annotate("text", x = 6.4, y = -1.05, label = "read-based route", colour = "#127a55", fontface = "bold", size = 3.6) +
        coord_cartesian(xlim = c(0.4, 8.6), ylim = c(-1.2, 1.2)) +
        theme_void()
    })

    seed_o <- reactiveVal(121)
    observeEvent(input$resample_o, seed_o(sample.int(1e6, 1)))
    shotgun <- reactive({
      set.seed(seed_o())
      ab <- c(input$abA, input$abB, input$abC); ab <- ab / sum(ab)
      w <- ab * mg_genomes$kb; w <- w / sum(w)                       # reads ∝ cells × genome size
      src <- sample(3, input$n_reads, replace = TRUE, prob = w)
      do.call(rbind, lapply(1:3, function(g) {
        G <- mg_genomes$kb[g] * 1000; n <- sum(src == g)
        if (n == 0) return(data.frame(g = g, start = numeric(0)))
        data.frame(g = g, start = sample.int(G - input$read_len + 1, n, replace = TRUE))
      }))
    })

    cov_profile <- reactive({
      r <- shotgun(); rl <- input$read_len
      do.call(rbind, lapply(1:3, function(g) {
        G <- mg_genomes$kb[g] * 1000
        d <- numeric(G + 1)
        s <- r$start[r$g == g]
        if (length(s)) { d <- d + tabulate(s, G + 1) - tabulate(pmin(s + rl, G + 1), G + 1) }
        data.frame(g = g, pos = seq_len(G), cov = cumsum(d)[seq_len(G)])
      }))
    })

    output$coverage_plot <- renderPlot({
      cp <- cov_profile()
      cp$name <- factor(mg_genomes$name[cp$g], levels = mg_genomes$name)
      # thin to ~600 points per genome for plotting
      cp <- do.call(rbind, lapply(split(cp, cp$g), function(d) d[unique(round(seq(1, nrow(d), length.out = 600))), ]))
      gaps <- cp[cp$cov == 0, ]
      means <- aggregate(cov ~ name, cp, mean)
      ggplot(cp, aes(pos / 1000, cov)) +
        geom_area(fill = "#cde2fb", colour = sv_col$primary, linewidth = 0.4) +
        geom_hline(data = means, aes(yintercept = cov), linetype = "22", colour = sv_col$ink2) +
        geom_rug(data = gaps, sides = "b", colour = sv_col$accent, alpha = 0.6, length = unit(0.06, "npc")) +
        facet_wrap(~name, ncol = 1, scales = "free_x") +
        labs(x = "Position (kb)", y = "Coverage") +
        theme_statviz() + theme(strip.text = element_text(face = "bold", colour = sv_col$ink, hjust = 0))
    })

    output$coverage_table <- renderUI({
      cp <- cov_profile(); r <- shotgun()
      ab <- c(input$abA, input$abB, input$abC); ab <- ab / sum(ab)
      rows <- lapply(1:3, function(g) {
        n <- sum(r$g == g); G <- mg_genomes$kb[g] * 1000; c <- n * input$read_len / G
        seen <- mean(cp$cov[cp$g == g] > 0)
        tags$tr(tags$td(mg_genomes$name[g]), tags$td(sprintf("%.0f%%", 100 * ab[g])),
                tags$td(sprintf("%.0f%%", 100 * n / nrow(r))), tags$td(fmt(c, 1)),
                tags$td(sprintf("%.0f%%", 100 * seen)), tags$td(sprintf("%.0f%%", 100 * (1 - exp(-c)))))
      })
      tagList(
        tags$table(class = "table table-sm small",
                   tags$thead(tags$tr(tags$th("Species"), tags$th("Cells"), tags$th("Reads"), tags$th("Coverage"),
                                      tags$th("Seen"), tags$th("Predicted"))),
                   tags$tbody(rows)),
        p(class = "small", strong("Cells vs reads: "), "species A has the biggest genome, so it gets a larger share of",
          "reads than of cells. ", strong("Seen vs predicted: "), "the fraction of each genome covered by at least",
          "one read, compared with the 1 − e^(−c) formula.")
      )
    })

    # ---- 2. reads & quality -------------------------------------------------
    qual <- reactive({
      set.seed(7)
      n <- 300; L <- 150; pos <- seq_len(L)
      base <- 37 - 14 * (pos / L)^3
      bad <- runif(n) < 0.06
      Q <- t(vapply(seq_len(n), function(i) {
        q <- base - if (bad[i]) 14 else 0
        q <- q + rnorm(L, 0, 2 + 5 * (pos / L)^2)
        dips <- which(runif(L) < 0.004); for (d in dips) q[d:min(L, d + 3)] <- q[d:min(L, d + 3)] - 15
        pmin(41, pmax(2, round(q)))
      }, numeric(L)))
      Q
    })

    trimmed <- reactive({
      Q <- qual(); thr <- input$q_thresh; L <- ncol(Q)
      len <- apply(Q, 1, function(q) {
        w <- stats::filter(q, rep(1 / 4, 4), sides = 1)[4:L]   # mean of bases i-3..i
        bad <- which(w < thr)
        if (length(bad)) max(0, bad[1] + 3 - 4) else L
      })
      list(len = len, keep = len >= input$min_len)
    })

    output$qual_plot <- renderPlot({
      Q <- qual()
      qs <- apply(Q, 2, quantile, c(0.05, 0.25, 0.5, 0.75, 0.95))
      df <- data.frame(pos = seq_len(ncol(Q)), q05 = qs[1, ], q25 = qs[2, ], q50 = qs[3, ], q75 = qs[4, ], q95 = qs[5, ])
      ggplot(df, aes(pos)) +
        annotate("rect", xmin = -Inf, xmax = Inf, ymin = 30, ymax = Inf, fill = "#e6f6ef", alpha = 0.6) +
        geom_ribbon(aes(ymin = q05, ymax = q95), fill = "#cde2fb") +
        geom_ribbon(aes(ymin = q25, ymax = q75), fill = "#86b6ef") +
        geom_line(aes(y = q50), colour = sv_col$primary, linewidth = 1.1) +
        geom_hline(yintercept = input$q_thresh, colour = sv_col$accent, linewidth = 1) +
        annotate("text", x = 3, y = 39.5, label = "Q ≥ 30: excellent", hjust = 0, size = 3.4, colour = "#127a55") +
        scale_y_continuous(limits = c(0, 41)) +
        labs(x = "Position in read", y = "Quality (Q)") + theme_statviz()
    })

    output$qc_boxes <- renderUI({
      t <- trimmed(); L <- ncol(qual())
      bases_kept <- sum(t$len[t$keep]) / (length(t$len) * L)
      perr <- 10^(-input$q_thresh / 10)
      layout_column_wrap(
        width = 1 / 3, fill = FALSE,
        sv_box("Reads kept", sprintf("%.0f%%", 100 * mean(t$keep)), tone = "blue"),
        sv_box("Bases kept", sprintf("%.0f%%", 100 * bases_kept), tone = "blue"),
        sv_box(sprintf("Q%d means", input$q_thresh), sprintf("1 in %s", format(round(1 / perr), big.mark = ",")),
               p("bases wrong"), tone = "gray")
      )
    })

    output$trim_hist <- renderPlot({
      t <- trimmed()
      ggplot(data.frame(len = t$len, keep = t$keep), aes(len, fill = keep)) +
        geom_histogram(binwidth = 5, boundary = 0, colour = sv_col$surface) +
        geom_vline(xintercept = input$min_len, linetype = "22", colour = sv_col$ink) +
        scale_fill_manual(values = c(`TRUE` = sv_col$primary, `FALSE` = sv_col$accent), guide = "none") +
        labs(x = "Length after trimming", y = "Reads") + theme_statviz(11)
    })

    # ---- 3. alignment -------------------------------------------------------
    aln <- reactive({
      clean <- function(s) substr(gsub("[^A-Z]", "", toupper(s)), 1, 25)
      a <- clean(input$seq_a); b <- clean(input$seq_b)
      validate(need(nchar(a) >= 1 && nchar(b) >= 1, "Type two sequences (letters only)."))
      dp_align(a, b, input$s_match, input$s_mis, input$s_gap, input$align_type)
    })

    output$dp_plot <- renderPlot({
      r <- aln(); S <- r$S
      df <- expand.grid(i = seq_len(nrow(S)), j = seq_len(ncol(S)))
      df$s <- S[cbind(df$i, df$j)]
      path <- as.data.frame(r$path); names(path) <- c("i", "j")
      rows <- c("–", r$x); cols <- c("–", r$y)
      small <- max(nrow(S), ncol(S)) <= 18
      g <- ggplot(df, aes(j, i)) +
        geom_tile(aes(fill = s), colour = sv_col$surface, linewidth = 0.6) +
        geom_tile(data = path, fill = NA, colour = sv_col$accent, linewidth = 1.3)
      if (small) g <- g + geom_text(aes(label = s, colour = s > max(S) * 0.6), size = 3.3) +
        scale_colour_manual(values = c(`TRUE` = "#ffffff", `FALSE` = sv_col$ink), guide = "none")
      g + scale_fill_gradientn(colours = seq_blues, name = "score") +
        scale_x_continuous(breaks = seq_along(cols), labels = cols, position = "top", expand = c(0, 0)) +
        scale_y_reverse(breaks = seq_along(rows), labels = rows, expand = c(0, 0)) +
        coord_equal() + labs(x = NULL, y = NULL) +
        theme_statviz() + theme(panel.grid = element_blank(), axis.text = element_text(size = 12, face = "bold"))
    })

    output$aln_text <- renderUI({
      r <- aln()
      tagList(
        tags$pre(class = "fs-5 mb-2", style = "line-height: 1.35;",
                 paste(r$a, r$mid, r$b, sep = "\n")),
        p(class = "small text-muted mb-0", "| = match, . = mismatch, – = gap.",
          if (input$align_type == "local") "Only the best-matching stretch is shown (local alignment).")
      )
    })

    output$aln_boxes <- renderUI({
      r <- aln()
      layout_column_wrap(
        width = 1 / 2, fill = FALSE,
        sv_box("Alignment score", fmt(r$score, 0), tone = "blue"),
        sv_box("Identity", sprintf("%.0f%%", 100 * r$identity), p(sprintf("over %d columns", nchar(r$a))), tone = "green")
      )
    })

    # ---- 4. assembly --------------------------------------------------------
    seed_g <- reactiveVal(131)
    observeEvent(input$new_genome, seed_g(sample.int(1e6, 1)))
    genome <- reactive({
      set.seed(seed_g())
      L <- input$g_len; R <- input$rep_len
      if (R == 0 || 2 * R > L - 10) return(rand_dna(L))
      rep <- rand_dna(R); rest <- L - 2 * R
      a <- round(rest * 0.3); b <- round(rest * 0.4)
      paste0(rand_dna(a), rep, rand_dna(b), rep, rand_dna(rest - a - b))
    })
    reads <- reactive({
      validate(need(input$a_read_len < nchar(genome()), "Reads must be shorter than the genome."))
      set.seed(seed_g() + 1)
      sample_reads(genome(), input$a_read_len, input$a_cov, input$a_err / 100)
    })
    kcounts <- reactive({
      validate(need(input$k <= input$a_read_len, "k must be no longer than the reads."))
      count_kmers(reads()$read, input$k)
    })
    true_kmers <- reactive({
      g <- genome(); k <- input$k; L <- nchar(g)
      unique(substring(g, 1:(L - k + 1), k:L))
    })
    asm <- reactive({
      kc <- kcounts(); keep <- names(kc)[kc >= input$min_count]
      validate(need(length(keep) > 0, "No k-mers left: lower the cut-off or raise the coverage."))
      debruijn_assemble(keep)
    })

    output$pileup <- renderPlot({
      r <- reads(); g <- genome(); rl <- input$a_read_len
      r <- r[order(r$start), ]
      if (nrow(r) > 250) r <- r[sort(sample(nrow(r), 250)), ]
      row_end <- numeric(0); r$row <- 0
      for (i in seq_len(nrow(r))) {
        free <- which(row_end < r$start[i])
        if (length(free)) { r$row[i] <- free[1] } else { row_end <- c(row_end, 0); r$row[i] <- length(row_end) }
        row_end[r$row[i]] <- r$start[i] + rl
      }
      err <- do.call(rbind, lapply(seq_len(nrow(r)), function(i) {
        true <- substr(g, r$start[i], r$start[i] + rl - 1)
        d <- which(chars(true) != chars(r$read[i]))
        if (length(d)) data.frame(x = r$start[i] + d - 1, row = r$row[i]) else NULL
      }))
      p <- ggplot(r) +
        geom_segment(aes(x = start, xend = start + rl - 1, y = row, yend = row), colour = sv_col$primary,
                     linewidth = 1.6, alpha = 0.8)
      if (!is.null(err)) p <- p + geom_point(data = err, aes(x, row), colour = sv_col$accent, shape = 15, size = 1.8)
      p + scale_y_reverse() + coord_cartesian(xlim = c(1, nchar(g))) +
        labs(x = "Position in genome", y = NULL) +
        theme_statviz(11) + theme(axis.text.y = element_blank(), panel.grid.major.y = element_blank())
    })

    output$kmer_hist <- renderPlot({
      kc <- kcounts()
      df <- data.frame(count = as.integer(kc), real = names(kc) %in% true_kmers())
      df$real <- factor(ifelse(df$real, "In the genome", "Error k-mer"), levels = c("In the genome", "Error k-mer"))
      ggplot(df, aes(count, fill = real)) +
        geom_bar(width = 0.8) +
        geom_vline(xintercept = input$min_count - 0.5, linetype = "22", colour = sv_col$ink) +
        scale_fill_manual(values = c(`In the genome` = sv_col$primary, `Error k-mer` = sv_col$accent), name = NULL) +
        labs(x = "Times seen", y = "Distinct k-mers") + theme_statviz(11)
    })

    output$asm_graph <- renderPlot({
      a <- asm(); n <- length(a$contigs)
      validate(need(n <= 250, sprintf("%d contigs: too many to draw. Try a larger k, higher coverage or a k-mer cut-off.", n)))
      if (n == 1) {
        return(ggplot() +
          annotate("point", x = 0, y = 0, size = 30, colour = sv_col$primary, alpha = 0.85) +
          annotate("text", x = 0, y = 0, label = nchar(a$contigs), colour = "#ffffff", size = 6, fontface = "bold") +
          annotate("text", x = 0, y = -0.55, size = 4.5, colour = sv_col$ink,
                   label = "No forks: every k-mer chains into a single contig.\nTry adding a repeat or sequencing errors.") +
          coord_cartesian(xlim = c(-1, 1), ylim = c(-1, 0.6)) + theme_void())
      }
      g <- igraph::graph_from_data_frame(a$edges, directed = TRUE, vertices = data.frame(id = seq_len(n)))
      set.seed(1)
      lay <- if (n == 1) matrix(c(0, 0), 1) else igraph::layout_nicely(g)
      nodes <- data.frame(id = seq_len(n), x = lay[, 1], y = lay[, 2], len = nchar(a$contigs))
      e <- a$edges
      if (nrow(e)) { e$x <- nodes$x[e$from]; e$y <- nodes$y[e$from]; e$xend <- nodes$x[e$to]; e$yend <- nodes$y[e$to]
        e <- e[e$from != e$to, ] }
      p <- ggplot()
      if (nrow(e)) p <- p + geom_segment(data = e, aes(x, y, xend = xend, yend = yend), colour = sv_col$muted,
                                         arrow = arrow(length = unit(0.15, "cm"), type = "closed"), linewidth = 0.5)
      p + geom_point(data = nodes, aes(x, y, size = len), colour = sv_col$primary, alpha = 0.85) +
        geom_text(data = nodes[nodes$len >= stats::quantile(nodes$len, 0.6) | n <= 25, ], aes(x, y, label = len),
                  colour = sv_col$ink, size = 3.2, vjust = -1.3) +
        scale_size_area(max_size = 14, guide = "none") +
        scale_x_continuous(expand = expansion(mult = 0.15)) + scale_y_continuous(expand = expansion(mult = 0.15)) +
        theme_void()
    })

    output$asm_boxes <- renderUI({
      a <- asm(); len <- nchar(a$contigs); G <- nchar(genome())
      layout_column_wrap(
        width = 1 / 2, fill = FALSE,
        sv_box("Contigs", length(len), p(sprintf("from %d distinct k-mers", a$n_kmers)), tone = "blue"),
        sv_box("Longest contig", sprintf("%.0f%%", 100 * min(1, max(len) / G)),
               p(sprintf("of the genome (%d of %d bases)", max(len), G),
                 if (max(len) < G && max(len) >= 0.85 * G) span(class = "small", br(), "The last few bases at each end are",
                                                                 "covered by fewer reads and often missed.")),
               tone = if (max(len) >= 0.9 * G) "green" else "orange"),
        sv_box("N50", n50(len), p("bases"), tone = "gray"),
        sv_box("Forks in the graph", a$branching, p("places with more than one way forward"), tone = "gray")
      )
    })

    output$asm_map <- renderPlot({
      a <- asm(); g <- genome(); G <- nchar(g)
      ct <- a$contigs[order(-nchar(a$contigs))]
      pos <- vapply(ct, function(s) as.integer(regexpr(s, g, fixed = TRUE)), 0L)
      placed <- data.frame(start = pos, len = nchar(ct))[pos > 0, ]
      placed <- head(placed, 25)
      placed$row <- seq_len(nrow(placed))
      ggplot() +
        annotate("segment", x = 1, xend = G, y = 0, yend = 0, colour = sv_col$good, linewidth = 4) +
        geom_segment(data = placed, aes(x = start, xend = start + len - 1, y = row, yend = row),
                     colour = sv_col$primary, linewidth = 3) +
        scale_y_reverse() + labs(x = "Position in genome", y = NULL) +
        theme_statviz(11) + theme(axis.text.y = element_blank(), panel.grid.major.y = element_blank())
    })

    # ---- 5. binning ---------------------------------------------------------
    seed_b <- reactiveVal(141)
    observeEvent(input$resample_b, seed_b(sample.int(1e6, 1)))
    contigs <- reactive(sim_contigs(input$n_gen, input$similar, input$contig_len, seed_b()))
    bins <- reactive({
      d <- contigs()
      X <- scale(cbind(d$gc, log10(pmax(d$coverage, 0.1))))
      set.seed(1)
      km <- kmeans(X, centers = input$n_gen, nstart = 25)
      # name each bin after its dominant genome where possible, for readable colours
      dom <- tapply(d$genome, km$cluster, function(v) names(which.max(table(v))))
      ord <- order(match(dom, sort(unique(d$genome))))
      d$bin <- factor(paste("Bin", match(km$cluster, ord)))
      d
    })
    bin_cols <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100")
    bin_shapes <- c(16, 17, 15, 18)

    bin_scatter <- function(d, var, lab) {
      lv <- levels(factor(d[[var]]))
      ggplot(d, aes(gc, coverage, colour = .data[[var]], shape = .data[[var]], size = length)) +
        geom_point(alpha = 0.75) +
        scale_colour_manual(values = setNames(bin_cols[seq_along(lv)], lv), name = lab) +
        scale_shape_manual(values = setNames(bin_shapes[seq_along(lv)], lv), name = lab) +
        scale_size_area(max_size = 6, guide = "none") +
        scale_y_log10() + scale_x_continuous(labels = function(v) paste0(round(100 * v), "%")) +
        labs(x = "GC content", y = "Coverage") + theme_statviz()
    }
    output$bin_truth <- renderPlot(bin_scatter(transform(contigs(), genome = factor(genome)), "genome", NULL))
    output$bin_pred <- renderPlot(bin_scatter(bins(), "bin", NULL))

    output$bin_table <- renderUI({
      d <- bins()
      tot <- tapply(d$length, d$genome, sum)
      rows <- lapply(levels(d$bin), function(b) {
        db <- d[d$bin == b, ]
        by_g <- tapply(db$length, db$genome, sum); by_g[is.na(by_g)] <- 0
        dom <- names(which.max(by_g))
        comp <- by_g[[dom]] / tot[[dom]]; cont <- 1 - by_g[[dom]] / sum(db$length)
        qual <- if (comp > 0.9 && cont < 0.05) "High quality" else if (comp >= 0.5 && cont < 0.1) "Medium quality" else "Low quality"
        tags$tr(tags$td(b), tags$td(nrow(db)), tags$td(format(sum(db$length), big.mark = ",")), tags$td(dom),
                tags$td(sprintf("%.0f%%", 100 * comp)), tags$td(sprintf("%.0f%%", 100 * cont)), tags$td(strong(qual)))
      })
      tagList(
        tags$table(class = "table table-sm",
                   tags$thead(tags$tr(tags$th("Bin"), tags$th("Contigs"), tags$th("Bases"), tags$th("Mostly from"),
                                      tags$th("Completeness"), tags$th("Contamination"), tags$th("MAG quality"))),
                   tags$tbody(rows)),
        p(class = "small text-muted mb-0",
          "Quality labels follow the MIMAG standard: high = more than 90% complete and under 5% contamination;",
          "medium = at least 50% complete and under 10% contamination.")
      )
    })

    # ---- 6. annotation ------------------------------------------------------
    db <- local({
      set.seed(2024)
      fun <- c("Beta-lactamase (antibiotic resistance)", "Nitrogenase iron protein (nitrogen fixation)",
               "Glycoside hydrolase (breaks down plant fibre)", "Methyl-coenzyme M reductase (methane production)",
               "ABC transporter ATP-binding protein (moves nutrients into the cell)",
               "Flagellin (builds the flagellum for swimming)")
      genes <- vapply(c(110, 95, 130, 85, 120, 100), rand_gene, "")
      data.frame(fun = fun, gene = genes, prot = vapply(genes, function(g) sub("\\*$", "", translate(g)), ""),
                 stringsAsFactors = FALSE)
    })

    seed_n <- reactiveVal(151)
    observeEvent(input$resample_n, seed_n(sample.int(1e6, 1)))
    contig <- reactive({
      set.seed(seed_n())
      pick <- sample(nrow(db), 3)
      parts <- character(0); truth <- data.frame(start = integer(0), end = integer(0), fun = character(0))
      pos <- 1
      for (i in seq_along(pick)) {
        spacer <- rand_dna(sample(150:350, 1)); parts <- c(parts, spacer); pos <- pos + nchar(spacer)
        g <- mutate_dna(db$gene[pick[i]], input$diverge / 100)
        substr(g, 1, 3) <- "ATG"                                 # keep the start codon
        parts <- c(parts, g)
        truth <- rbind(truth, data.frame(start = pos, end = pos + nchar(g) - 1, fun = db$fun[pick[i]]))
        pos <- pos + nchar(g)
      }
      s <- paste(parts, collapse = "")
      s <- paste0(s, rand_dna(max(0, 2000 - nchar(s))))
      list(seq = substr(s, 1, 2000), truth = truth[truth$end <= 2000, ])
    })

    # ungapped "seed and extend"-style comparison: best diagonal identity
    best_hit <- function(p) {
      x <- chars(p)
      res <- t(vapply(db$prot, function(q) {
        y <- chars(q); best <- c(0, 0)
        for (off in seq(-(length(x) - 20), length(y) - 20)) {
          i <- seq_along(x); j <- i + off; ok <- j >= 1 & j <= length(y)
          if (sum(ok) < 20) next
          m <- sum(x[ok] == y[j[ok]])
          if (m > best[1]) best <- c(m, sum(ok))
        }
        c(best[1], best[2])
      }, numeric(2)))
      k <- which.max(res[, 1])
      list(fun = db$fun[k], ident = res[k, 1] / max(1, res[k, 2]), len = res[k, 2])
    }

    orfs <- reactive({
      c0 <- contig(); o <- find_orfs(c0$seq, input$min_orf)
      if (!nrow(o)) return(o)
      o <- o[order(o$start), ]
      o$prot <- vapply(seq_len(nrow(o)), function(i) sub("\\*$", "", translate(substr(c0$seq, o$start[i], o$end[i]))), "")
      hits <- lapply(o$prot, best_hit)
      o$ident <- vapply(hits, `[[`, 0, "ident"); o$alen <- vapply(hits, `[[`, 0, "len")
      o$hit <- o$ident >= 0.3 & o$alen >= 30
      o$fun <- ifelse(o$hit, vapply(hits, `[[`, "", "fun"), "Hypothetical protein (no match)")
      o$id <- paste0("ORF", seq_len(nrow(o)))
      o
    })

    output$orf_pick <- renderUI({
      o <- orfs()
      if (nrow(o)) selectInput(session$ns("orf_sel"), "Show the codons of", o$id)
    })

    output$orf_map <- renderPlot({
      o <- orfs(); tr <- contig()$truth
      g <- ggplot() +
        geom_rect(data = tr, aes(xmin = start, xmax = end, ymin = -0.25, ymax = 0.25), fill = "#bfe9d6") +
        annotate("text", x = 5, y = 0, label = "real genes", hjust = 0, size = 3.2, colour = "#127a55")
      if (nrow(o)) {
        g <- g + geom_segment(data = o, aes(x = start, xend = end, y = frame, yend = frame, colour = hit),
                              linewidth = 5, arrow = arrow(length = unit(0.25, "cm"), type = "closed")) +
          geom_text(data = o, aes(x = (start + end) / 2, y = frame + 0.38, label = id), size = 3.2, colour = sv_col$ink)
      }
      g + scale_colour_manual(values = c(`TRUE` = sv_col$primary, `FALSE` = sv_col$muted), guide = "none") +
        scale_y_continuous(breaks = 0:3, labels = c("", "Frame 1", "Frame 2", "Frame 3"), limits = c(-0.4, 3.6)) +
        coord_cartesian(xlim = c(1, 2000)) +
        labs(x = "Position on contig", y = NULL) + theme_statviz(11)
    })

    output$annot_table <- renderUI({
      o <- orfs()
      if (!nrow(o)) return(p("No open reading frames this long. Lower the minimum gene length."))
      tagList(
        tags$table(class = "table table-sm small",
                   tags$thead(tags$tr(tags$th("Gene"), tags$th("Position"), tags$th("Length (aa)"),
                                      tags$th("Best match"), tags$th("Identity"))),
                   tags$tbody(lapply(seq_len(nrow(o)), function(i) tags$tr(
                     tags$td(o$id[i]), tags$td(sprintf("%d-%d", o$start[i], o$end[i])), tags$td(nchar(o$prot[i])),
                     tags$td(if (o$hit[i]) strong(o$fun[i]) else span(class = "text-muted", o$fun[i])),
                     tags$td(sprintf("%.0f%%", 100 * o$ident[i])))))),
        p(class = "small text-muted mb-0",
          "Database: six made-up proteins named after real functions. A match needs at least 30% identity over",
          "30 or more amino acids. Spurious ORFs in non-coding DNA match nothing.")
      )
    })

    output$codon_view <- renderUI({
      o <- orfs()
      req(nrow(o), input$orf_sel %in% o$id)
      i <- match(input$orf_sel, o$id)
      dna <- substr(contig()$seq, o$start[i], o$end[i])
      n <- min(12, nchar(dna) %/% 3)
      cod <- substring(dna, seq(1, by = 3, length.out = n), seq(3, by = 3, length.out = n))
      tagList(
        p(class = "small", "The first codons of ", strong(input$orf_sel), ". Each three DNA letters code for one",
          "amino acid (shown by its one-letter code)."),
        tags$table(class = "table table-sm text-center font-monospace",
                   tags$tr(lapply(cod, function(c) tags$td(c))),
                   tags$tr(lapply(cod, function(c) tags$td(strong(CODON_TABLE[[c]]))))),
        p(class = "small", "Protein (", nchar(o$prot[i]), " amino acids):"),
        tags$pre(class = "small", style = "white-space: pre-wrap; word-break: break-all;", o$prot[i])
      )
    })

    # ---- 7. taxonomy --------------------------------------------------------
    output$tax_sliders <- renderUI({
      ns <- session$ns; def <- c(40, 25, 15, 10, 10)
      tagList(lapply(seq_len(nrow(taxa)), function(i)
        sliderInput(ns(paste0("tax", i)), taxa$taxon[i], 0, 100, def[i], step = 1)))
    })

    profile <- reactive({
      ab <- vapply(seq_len(nrow(taxa)), function(i) input[[paste0("tax", i)]] %||% c(40, 25, 15, 10, 10)[i], 0)
      validate(need(sum(ab) > 0, "Give at least one taxon a non-zero share."))
      ab <- ab / sum(ab)
      set.seed(5)
      share_reads <- ab * taxa$genome_mb / sum(ab * taxa$genome_mb)
      reads <- as.vector(rmultinom(1, input$t_reads, share_reads))
      n <- nrow(taxa)
      assigned <- numeric(n); unclassified <- 0
      for (i in seq_len(n)) {
        if (taxa$in_db[i]) {
          a <- rbinom(1, reads[i], 0.93); assigned[i] <- assigned[i] + a; unclassified <- unclassified + reads[i] - a
        } else {
          wrong <- rbinom(1, reads[i], 0.15); assigned[4] <- assigned[4] + wrong; unclassified <- unclassified + reads[i] - wrong
        }
      }
      est <- assigned / taxa$genome_mb; est <- est / sum(est)
      list(ab = ab, reads = reads, share_reads = reads / sum(reads), assigned = assigned, est = est,
           unclassified = unclassified / sum(reads))
    })

    output$tax_bars <- renderPlot({
      pr <- profile()
      short <- c("B. fragilis", "E. coli", "L. reuteri", "P. copri", "Unknown")
      df <- rbind(data.frame(taxon = short, what = "True share of cells", v = pr$ab),
                  data.frame(taxon = short, what = "Share of reads", v = pr$share_reads),
                  data.frame(taxon = short, what = "Estimate (size-corrected)", v = pr$est))
      df$what <- factor(df$what, levels = c("True share of cells", "Share of reads", "Estimate (size-corrected)"))
      df$taxon <- factor(df$taxon, levels = short)
      ggplot(df, aes(taxon, v, fill = what)) +
        geom_col(position = position_dodge(0.8), width = 0.75) +
        scale_fill_manual(values = c(`True share of cells` = sv_col$good, `Share of reads` = sv_col$accent,
                                     `Estimate (size-corrected)` = sv_col$primary), name = NULL) +
        scale_y_continuous(labels = function(v) paste0(100 * v, "%"), expand = expansion(mult = c(0, 0.05))) +
        labs(x = NULL, y = "Share") + theme_statviz()
    })

    output$tax_boxes <- renderUI({
      pr <- profile()
      H_true <- -sum(ifelse(pr$ab > 0, pr$ab * log(pr$ab), 0))
      H_est <- -sum(ifelse(pr$est > 0, pr$est * log(pr$est), 0))
      tagList(
        sv_box("Reads unclassified", sprintf("%.0f%%", 100 * pr$unclassified),
               p("mostly from the species missing from the database"), tone = "gray"),
        sv_box("Shannon diversity", fmt(H_est), p(sprintf("estimated (true value %s; maximum for 5 taxa: %s)",
                                                           fmt(H_true), fmt(log(5)))), tone = "blue")
      )
    })

    kmer_case <- reactive({
      set.seed(3)
      n <- 30
      lab <- switch(
        input$kmer_case,
        clear = sample(c("Escherichia coli", "Escherichia (genus, shared)", "No match"), n, TRUE, c(0.75, 0.2, 0.05)),
        shared = sample(c("Escherichia coli", "Escherichia (genus, shared)", "No match"), n, TRUE, c(0.1, 0.8, 0.1)),
        novel = sample(c("Escherichia coli", "Escherichia (genus, shared)", "No match"), n, TRUE, c(0.03, 0.05, 0.92))
      )
      data.frame(i = seq_len(n), lab = factor(lab, levels = c("Escherichia coli", "Escherichia (genus, shared)", "No match")))
    })

    output$kmer_read <- renderPlot({
      d <- kmer_case()
      ggplot(d, aes(i, 1, fill = lab)) +
        geom_tile(width = 0.9, height = 0.7, show.legend = TRUE) +
        scale_fill_manual(values = c(`Escherichia coli` = sv_col$primary, `Escherichia (genus, shared)` = sv_col$good,
                                     `No match` = "#d6d5d0"), name = "k-mer found in", drop = FALSE) +
        labs(x = NULL, y = NULL) + theme_void() + theme(legend.position = "right")
    })

    output$kmer_verdict <- renderUI({
      d <- kmer_case(); t <- table(d$lab); n <- nrow(d)
      verdict <- switch(
        input$kmer_case,
        clear = "Many k-mers are unique to E. coli, so the read is assigned to the species Escherichia coli.",
        shared = "Most matching k-mers are shared by several Escherichia species, so the read is assigned to the genus Escherichia (their lowest common ancestor): correct, just less specific.",
        novel = "Almost no k-mers match the database, so the read is left unclassified. The organism may be new to science or simply missing from the database."
      )
      p(class = "small mx-3 mt-2", sprintf("%d of %d k-mers unique to E. coli, %d shared within the genus, %d with no match. ",
                                          t[[1]], n, t[[2]], t[[3]]), strong(verdict))
    })

    quiz_server(input, output, session, list(
      list(q = "Shotgun metagenomics sequences...",
           choices = c("one gene from every organism", "random fragments of all the DNA in a sample",
                       "only organisms that can be grown in the lab", "the RNA of a single species"),
           answer = 2, why = "DNA from every cell is fragmented and sequenced at random, giving reads from all genomes mixed together."),
      list(q = "A base with quality Q30 has an error probability of about...",
           choices = c("1 in 30", "1 in 100", "1 in 1,000", "30%"),
           answer = 3, why = "p = 10^(−Q/10) = 10^(−3) = 0.001."),
      list(q = "In a de Bruijn graph assembler, what usually breaks a contig?",
           choices = c("A repeat longer than k", "A long read", "High coverage", "A high GC content"),
           answer = 1, why = "Repeats longer than k create forks: the graph can't tell which copy leads where."),
      list(q = "Why do sequencing errors show up as k-mers seen only once or twice?",
           choices = c("Errors are always at the end of reads",
                       "An error is random, so the false k-mers it creates are unlikely to be repeated in other reads",
                       "The sequencer flags them", "They aren't; errors look like real k-mers"),
           answer = 2, why = "Real k-mers are seen about 'coverage' times; error k-mers are mostly unique. Ignoring rare k-mers removes most errors."),
      list(q = "Binning groups contigs into genomes using mainly...",
           choices = c("contig length", "GC content / sequence composition and coverage",
                       "the order contigs were assembled in", "gene names"),
           answer = 2, why = "Pieces of one genome share its composition 'accent' and its abundance (coverage)."),
      list(q = "Species X has a genome twice as long as species Y, and they have the same number of cells. X will contribute...",
           choices = c("the same number of reads", "about twice as many reads", "half as many reads", "no reads"),
           answer = 2, why = "Reads are sampled in proportion to DNA, so read share must be divided by genome size to estimate cell share."),
      list(q = "Dynamic programming in alignment is useful because...",
           choices = c("it tries every possible alignment one by one",
                       "it builds the best alignment from the best solutions to smaller sub-problems, filling a table",
                       "it only works for identical sequences", "it doesn't allow gaps"),
           answer = 2, why = "Each cell needs only its three neighbours, so the optimal alignment is found without trying every possibility.")
    ))
  })
}
