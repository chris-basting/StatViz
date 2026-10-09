# Principal Component Analysis module ------------------------------------------
#
# Tabs:
#   1. The big idea   - rotate a line through a 2-D cloud; find max variance
#   2. Step by step   - center/scale -> covariance -> eigen -> scores
#   3. How many PCs?  - scree plot & cumulative variance
#   4. Interpreting   - loadings & biplot
#   5. Reconstruction - rebuild data from k components
#   6. Into 3D        - a 3-D cloud flattened onto its best plane
#   7. Quiz

data_tabs <- c("steps", "scree", "interp", "recon")

pca_ui <- function(id) {
  ns <- NS(id)
  is_tab <- function(tabs) {
    sprintf("[%s].includes(input.tab)", paste0("'", tabs, "'", collapse = ","))
  }

  layout_sidebar(
    sidebar = sidebar(
      width = 300,
      # ---- 2-D intuition controls
      conditionalPanel(
        is_tab("idea"), ns = ns,
        h6("Shape of the cloud"),
        sliderInput(ns("n2"), "Number of points", 20, 400, 120, step = 10),
        sliderInput(ns("rho"), "Correlation between x1 and x2", -0.95, 0.95, 0.75, step = 0.05),
        sliderInput(ns("sdx"), "Spread of x1 (SD)", 0.5, 3, 1.5, step = 0.1),
        sliderInput(ns("sdy"), "Spread of x2 (SD)", 0.5, 3, 1, step = 0.1),
        actionButton(ns("resample2"), "Draw a new sample", icon = icon("shuffle"),
                     class = "btn-outline-secondary btn-sm w-100"),
        hr(),
        h6("Your direction"),
        sliderInput(ns("angle"), "Angle θ of the line (degrees)", 0, 180, 20, step = 1,
                    post = "°"),
        helpText("θ is measured counter-clockwise from the x1 axis."),
        actionButton(ns("play"), "Rotate", icon = icon("play"),
                     class = "btn-outline-primary btn-sm w-100 mb-2"),
        tags$script(src = "rotate.js"),
        tags$script(HTML(sprintf(
          "svRotate({slider: '%s', button: '%s', plot: '%s', stop: ['%s']});",
          ns("angle"), ns("play"), ns("idea_scatter"), ns("snap")))),
        checkboxInput(ns("show_proj"), "Show projections onto the line", TRUE),
        checkboxInput(ns("show_resid"), "Show residuals (distances to the line)", TRUE),
        actionButton(ns("snap"), "Reveal PC1", icon = icon("bullseye"),
                     class = "btn-primary w-100")
      ),
      # ---- real-data controls
      conditionalPanel(
        is_tab(data_tabs), ns = ns,
        selectInput(ns("dataset"), "Dataset", pca_dataset_choices),
        conditionalPanel(
          "input.dataset == 'upload'", ns = ns,
          fileInput(ns("file"), NULL, accept = c(".csv", "text/csv"),
                    placeholder = "CSV with a header row")
        ),
        uiOutput(ns("blurb")),
        selectizeInput(ns("vars"), "Variables to include", choices = NULL, multiple = TRUE,
                       options = list(plugins = list("remove_button"))),
        selectInput(ns("group"), "Color points by", choices = c("(none)" = "")),
        checkboxInput(ns("scale"), "Scale variables to unit variance", TRUE),
        helpText("Scaling = PCA on the correlation matrix. Without it, variables",
                 "with big units dominate."),
        conditionalPanel(
          is_tab(c("steps", "interp")), ns = ns,
          hr(),
          layout_columns(
            selectInput(ns("pcx"), "X axis", "PC1"),
            selectInput(ns("pcy"), "Y axis", "PC2")
          )
        )
      ),
      # ---- 3-D controls
      conditionalPanel(
        is_tab("3d"), ns = ns,
        sliderInput(ns("n3"), "Number of points", 50, 500, 200, step = 25),
        sliderInput(ns("thick"), "Thickness of the cloud (SD off the plane)",
                    0.1, 3, 0.4, step = 0.1),
        actionButton(ns("resample3"), "New random tilt", icon = icon("shuffle"),
                     class = "btn-outline-secondary btn-sm w-100"),
        hr(),
        checkboxInput(ns("show_plane"), "Show the PC1-PC2 plane", TRUE),
        checkboxInput(ns("show_axes"), "Show PC arrows", TRUE)
      ),
      conditionalPanel(
        is_tab("quiz"), ns = ns,
        p("Test yourself on the ideas from the other tabs. Answers are checked",
          "when you press the button at the bottom.")
      )
    ),

    navset_underline(
      id = ns("tab"),

      # ==== 1. Big idea ========================================================
      nav_panel(
        "1 · The big idea", value = "idea",
        explain(withMathJax(
          p(strong("The goal:"), "describe each point with ", em("one"), " number instead of two,",
            "while losing as little information as possible. PCA does this by finding the single",
            "line through the data that captures the most spread."),
          p(strong("Try it:"), "drag the angle slider (or press ", strong("Rotate"), ").",
            "Each blue point drops straight onto the orange line, like casting a shadow.",
            "Where it lands (orange dot) is the point's ", strong("score"),
            ": one number saying where the point sits along the line.",
            "The gray segment is the point's ", strong("residual"),
            ": the part of the point that the line misses."),
          p(strong("Keeping track of what's kept and lost."),
            "To compare lines, square each score and each residual and add them up:"),
          math_words(
            "$$\\sum_i d_i^2 \\;=\\; \\sum_i s_i^2 \\;+\\; \\sum_i r_i^2$$",
            "For every point, the score and residual form a right triangle with the point's distance",
            "from the center as the long side, so Pythagoras says",
            "(distance)² = (score)² + (residual)². Adding that up over all points gives",
            strong("total spread = spread kept + spread lost."),
            "Each piece is called a ", strong("sum of squares (SS)"), ".",
            symbols = list(
              "\\(\\sum_i\\)" = "add up over every point",
              "\\(d_i\\)" = "distance of point i from the center",
              "\\(s_i\\)" = "its score",
              "\\(r_i\\)" = "its residual"
            )
          ),
          p(strong("Why that matters:"), "the total on the left is the same for every line, because",
            "the points never move. So the line that keeps the most automatically loses the least.",
            "That line is the ", strong("first principal component (PC1)"), ".",
            "You'll often hear it called 'the direction of maximum ", strong("variance"), "'.",
            "Variance is just the average squared distance from the center",
            "(the sum of squares divided by n − 1), so it picks the same winner."),
          p(class = "small text-muted mb-2",
            "Not the same as regression: a regression line measures its misses ", em("vertically"),
            " because it's trying to predict x2 from x1. PCA measures them ", em("perpendicular"),
            " to the line, treating both variables as equals.")
        )),
        layout_columns(
          col_widths = c(7, 5),
          card(
            full_screen = TRUE,
            card_header("The data and your line"),
            plotOutput(ns("idea_scatter"), height = "460px"),
            axes_note(
              "x1, the first of two variables in a simulated dataset (think height).",
              "x2, the second variable (think weight).",
              "Both are centered: each variable's average is subtracted, so the cloud sits at (0, 0)",
              "and every candidate line passes through the middle of the data. The orange arc marks",
              "the angle θ of your line; gray segments are the residuals."
            )
          ),
          div(
            uiOutput(ns("idea_stats")),
            card(
              card_header("Sums of squares vs. angle"),
              plotOutput(ns("idea_curve"), height = "280px"),
              axes_note(
                "the angle θ of the orange line, from 0° to 180° (every possible line).",
                paste("sum of squares. Blue = spread kept (squared scores, added up).",
                      "Gray = spread lost (squared residuals, added up)."),
                "Blue + gray always equals the dotted total. The blue peak and the gray dip happen at the",
                "same angle: PC1. PC2 sits exactly 90° away, where the line keeps the least."
              )
            )
          )
        )
      ),

      # ==== 2. Step by step ====================================================
      nav_panel(
        "2 · Step by step", value = "steps",
        explain(
          p("Turning a line by hand works for two variables, but real datasets have many.",
            "These five steps are the recipe a computer follows to find all the best directions at once.",
            "Each step shows the formula, what it means in plain words, and the result on real data",
            "(pick a dataset in the sidebar).")
        ),
        accordion(
          id = ns("steps_acc"), open = c("Step 1"), multiple = TRUE,
          accordion_panel(
            value = "Step 1", title = "Step 1 · Put every variable on the same footing",
            icon = icon("1"),
            withMathJax(
              math_words(
                "$$z_{ij} = \\frac{x_{ij} - \\bar{x}_j}{s_j}$$",
                "take each value, subtract its column's average, and (if scaling is on) divide by the",
                "column's standard deviation, which is its typical distance from the average.",
                "Afterwards every variable averages 0, so the cloud is centered on the origin, just like",
                "in tab 1. With scaling, every variable also has a spread of 1, so a variable measured",
                "in big units (dollars, say) can't dominate one measured in small units (percentages)",
                "just because its numbers are bigger.",
                symbols = list(
                  "\\(x_{ij}\\)" = "original value of variable j for observation i",
                  "\\(\\bar{x}_j\\)" = "average of variable j",
                  "\\(s_j\\)" = "standard deviation of variable j",
                  "\\(z_{ij}\\)" = "the new, centered (and scaled) value"
                )
              )
            ),
            layout_columns(
              col_widths = c(5, 7),
              div(p(class = "small text-muted", "Each variable before and after:"),
                  tableOutput(ns("step1_summary"))),
              div(p(class = "small text-muted", "First rows of the transformed data (called Z):"),
                  tableOutput(ns("step1_head")))
            )
          ),
          accordion_panel(
            value = "Step 2", title = "Step 2 · Measure which variables move together",
            icon = icon("2"),
            withMathJax(
              math_words(
                "$$S = \\frac{1}{n-1} Z^\\top Z$$",
                "for every pair of variables, multiply their centered values observation by observation,",
                "then average those products. If two variables tend to be above average at the same",
                "time, most products are positive, giving a positive ", strong("covariance"), ".",
                "If one tends to be high when the other is low, it's negative. The results fill a table",
                "S, the ", strong("covariance matrix"), ", with one row and one column per variable.",
                symbols = list(
                  "\\(Z\\)" = "the centered data table from step 1",
                  "\\(Z^\\top Z\\)" = "shorthand for 'multiply every pair of columns and add up'",
                  "\\(n\\)" = "number of observations"
                )
              ),
              p(strong("Why we care:"), "the diagonal of S holds each variable's own variance (its spread).",
                "With scaling these are all 1, and S becomes the familiar ", strong("correlation matrix"), ".",
                "Big numbers off the diagonal mean ", strong("redundancy"),
                ": two variables telling much the same story.",
                "That redundancy is what lets PCA summarize many variables with just a few.")
            ),
            plotOutput(ns("step2_cov"), height = "380px"),
            p(class = "small sv-axes", strong("Rows and columns: "),
              "the original variables. Each cell holds the covariance (or correlation, when scaled)",
              "of its row variable with its column variable. Red = move together, blue = move in",
              "opposite directions, gray = unrelated.")
          ),
          accordion_panel(
            value = "Step 3", title = "Step 3 · Why the covariance matrix holds the answer",
            icon = icon("3"),
            withMathJax(
              p("This is the key link. In tab 1 we found the best line by trying every angle.",
                "We don't have to do that, for two reasons."),
              p(strong("(a) S can tell us the spread in any direction.")),
              math_words(
                "$$\\text{spread along } u \\;=\\; u^\\top S\\, u$$",
                "pick a direction u: an arrow of length 1, like the orange line in tab 1.",
                "'Sandwiching' S between two copies of u gives the variance of the scores along that",
                "direction. That's exactly the blue curve in tab 1, but worked out from the small table S",
                "instead of from all the data points.",
                symbols = list(
                  "\\(u\\)" = "a direction (a pair of numbers for 2 variables, a list of p numbers for p variables)",
                  "\\(u^\\top S\\, u\\)" = "multiply u by S, then by u again; the result is one number"
                )
              ),
              deeper(
                "Show where this comes from",
                p("The scores along u are \\(Zu\\) (each observation's values weighted by u and added up).",
                  "Their variance is their sum of squares divided by \\(n-1\\):",
                  "$$\\operatorname{Var}(Zu) = \\frac{(Zu)^\\top (Zu)}{n-1} = u^\\top \\Big(\\frac{Z^\\top Z}{n-1}\\Big)\\, u = u^\\top S\\, u$$")
              ),
              p(strong("(b) The best direction is one that S doesn't turn."),
                "Picture S as a machine: feed it a direction u and it outputs a new arrow Su.",
                "For most directions the output points somewhere else, because S bends arrows toward where the",
                "data are most spread out. But a few special directions come out pointing the ",
                em("same way"), ", only longer or shorter. These are the ", strong("eigenvectors"),
                " of S ('eigen' is German for 'own': the matrix's own directions)."),
              math_words(
                "$$S\\,u = \\lambda\\,u$$",
                "putting u through S gives back the same direction, just stretched by some number λ",
                "(the Greek letter lambda). u is an ", strong("eigenvector"), ", and the stretch factor λ is its ",
                strong("eigenvalue"), ".",
                symbols = list("\\(\\lambda\\)" = "the stretch factor (a single number)")
              ),
              p(strong("The punchline:"), "the direction of maximum spread is always an eigenvector,",
                "and its eigenvalue λ is exactly the spread (variance) along it.",
                "So the eigenvector with the biggest eigenvalue is PC1, the next biggest is PC2, and so on."),
              deeper(
                "Show the calculus (why the best direction must be an eigenvector)",
                p("We want the largest \\(u^\\top S u\\) among arrows of length 1 (\\(u^\\top u = 1\\)).",
                  "A standard calculus tool for 'maximize subject to a constraint' is a Lagrange multiplier \\(\\lambda\\).",
                  "Setting the slope to zero:",
                  "$$\\frac{\\partial}{\\partial u}\\Big[u^\\top S u - \\lambda(u^\\top u - 1)\\Big] = 2Su - 2\\lambda u = 0",
                  "\\quad\\Longrightarrow\\quad S\\,u = \\lambda\\,u$$",
                  "And the spread along such a u is \\(u^\\top S u = u^\\top(\\lambda u) = \\lambda\\, u^\\top u = \\lambda\\).")
              ),
              p(strong("(c) Try it."), "Below is the 2-D cloud from tab 1. Turn the orange arrow u and watch",
                "the blue arrow Su. When they line up, you've found an eigenvector.")
            ),
            layout_columns(
              col_widths = c(7, 5),
              div(
                plotOutput(ns("bridge_plot"), height = "420px"),
                axes_note("x1 (centered), the first variable of the tab-1 cloud.",
                          "x2 (centered), the second variable.",
                          "Gray dots: the data. Gray circle: every possible direction u. Blue ellipse: where S sends",
                          "each of them. Orange arrow: u. Blue arrow: Su. Green dashed: the eigenvectors.")
              ),
              div(
                sliderInput(ns("u_angle"), "Direction u (angle)", 0, 180, 110, step = 1, post = "°"),
                actionButton(ns("power_step"), "Feed Su back in as the new u", icon = icon("forward-step"),
                             class = "btn-outline-primary btn-sm w-100"),
                p(class = "small text-muted mt-1",
                  "Click a few times. Each click takes the output Su and uses it as the next input.",
                  "Because S keeps bending arrows toward the direction of most spread, this homes in on PC1.",
                  "This trick, called ", em("power iteration"), ", is one simple way computers find eigenvectors."),
                uiOutput(ns("bridge_math"))
              )
            )
          ),
          accordion_panel(
            value = "Step 4", title = "Step 4 · Find all the eigenvectors at once",
            icon = icon("4"),
            withMathJax(
              p("With p variables, S has p eigenvectors, all at right angles to each other.",
                "Software finds them all in one go (R's ", code("eigen()"), " function).",
                "This is called an ", strong("eigen-decomposition"), "."),
              math_words(
                "$$S\\,v_k = \\lambda_k\\,v_k \\qquad k = 1, 2, \\dots, p$$",
                "each principal component k has a direction \\(v_k\\) and a spread \\(\\lambda_k\\).",
                "Sorted from biggest spread to smallest they're PC1, PC2, PC3, ...",
                "Because the directions are at right angles, each PC captures a different pattern with",
                "no overlap (the PCs are uncorrelated).",
                symbols = list(
                  "\\(v_k\\)" = "direction of PC k (its loadings)",
                  "\\(\\lambda_k\\)" = "variance captured by PC k",
                  "\\(p\\)" = "number of variables"
                )
              )
            ),
            layout_columns(
              col_widths = c(4, 8),
              div(
                p(class = "small text-muted",
                  "Each PC's eigenvalue (the spread it captures) and its share of the total.",
                  "'Var. of scores' measures that spread directly from the data, and it matches."),
                tableOutput(ns("step3_eigs"))
              ),
              plotOutput(ns("step3_vecs"), height = "360px")
            ),
            p(class = "small text-muted",
              "Right: the eigenvectors. Each column is a PC; each number (a ", strong("loading"), ")",
              "says how much an original variable contributes to it. The signs are arbitrary:",
              "flipping every sign in a column describes the same line, pointed the other way."),
            uiOutput(ns("step3_check")),
            p(strong("See for yourself:"), "multiplying S by each eigenvector gives the same numbers as",
              "multiplying the eigenvector by its λ. S only stretches these directions."),
            tableOutput(ns("step3_verify"))
          ),
          accordion_panel(
            value = "Step 5", title = "Step 5 · Describe every observation with the new axes",
            icon = icon("5"),
            withMathJax(
              math_words(
                "$$T = Z\\,V$$",
                "each observation's score on a PC is a weighted sum of its centered values, using",
                "that PC's loadings as the weights. Do this for every observation and every PC.",
                symbols = list(
                  "\\(Z\\)" = "centered data (step 1)",
                  "\\(V\\)" = "eigenvectors, one per column (step 4)",
                  "\\(T\\)" = "the scores: one row per observation, one column per PC"
                )
              )
            ),
            uiOutput(ns("step5_example")),
            p("Plotting the first two columns of scores gives the classic PCA plot: a 2-D map of the",
              "data that keeps as much of the spread as any 2-D view can.",
              "Hover over points to see which observation is which."),
            plotlyOutput(ns("step4_scores"), height = "460px"),
            div(class = "small sv-axes",
                div(strong("Horizontal / vertical axis: "),
                    "each observation's score on the chosen PC: where it sits along that new axis.",
                    "0 = an average observation."),
                div(class = "text-muted",
                    "The % in each axis title is the share of the total spread that PC captures.",
                    "Choose which PCs to plot in the sidebar."))
          )
        )
      ),

      # ==== 3. Scree ===========================================================
      nav_panel(
        "3 · How many PCs?", value = "scree",
        explain(
          p("There are as many PCs as there are variables. Each one captures a share of the total",
            "spread: PC1 the most, then PC2, and so on. Often the first few capture nearly everything,",
            "and the rest is mostly noise. Keeping only those few is how PCA simplifies data."),
          p("How many to keep is a judgment call. Three common rules of thumb:"),
          tags$ul(
            tags$li(strong("Enough to reach a target:"), " keep PCs until together they capture, say,",
                    "80-90% of the spread (the line in the chart)."),
            tags$li(strong("The elbow:"), " look for where the bars suddenly level off, and keep",
                    "the PCs before that point. 'Scree' is the rubble at the foot of a cliff:",
                    "the PCs after the elbow."),
            tags$li(strong("Kaiser rule"), " (scaled data only): keep PCs whose eigenvalue is above 1.",
                    "With scaling, every original variable has a spread of 1, so a PC below 1",
                    "summarizes less than a single original variable would on its own.")
          )
        ),
        layout_columns(
          col_widths = c(8, 4),
          card(
            card_header("Scree plot"),
            plotOutput(ns("scree_plot"), height = "400px"),
            axes_note("the principal components, in order of how much variance they explain.",
                      "% of the total variance. Bars = each PC on its own; the line = running total.")
          ),
          div(
            card(
              card_header("Your threshold"),
              sliderInput(ns("thresh"), "Variance to keep (%)", 50, 99, 80, step = 1),
              uiOutput(ns("scree_answer"))
            ),
            card(card_header("Eigenvalue table"), tableOutput(ns("scree_table")))
          )
        )
      ),

      # ==== 4. Interpreting ====================================================
      nav_panel(
        "4 · Interpreting", value = "interp",
        explain(
          p("Each PC is a recipe: a weighted mix of the original variables. The weights, called",
            strong("loadings"), ", tell you what the PC 'means'.",
            "Variables with large loadings of the same sign rise and fall together along that PC;",
            "opposite signs pull in opposite directions."),
          p("The ", strong("biplot"), " overlays both: points are observations (scores),",
            "arrows are variables (loadings). Arrows pointing the same way are positively correlated;",
            "at right angles, roughly uncorrelated; opposite, negatively correlated.")
        ),
        layout_columns(
          col_widths = c(7, 5),
          card(full_screen = TRUE, card_header("Biplot"),
               plotOutput(ns("biplot"), height = "520px"),
               axes_note("scores on the PC chosen as 'X axis' in the sidebar.",
                         "scores on the PC chosen as 'Y axis'.",
                         "Dots are observations. Arrows are variables: each arrow's tip is at",
                         "(loading on X PC, loading on Y PC), stretched by a common factor so they",
                         "are visible. Longer arrows = variables better represented by these two PCs.")),
          card(card_header("Loadings"),
               selectInput(ns("load_pc"), NULL, "PC1", width = "140px"),
               plotOutput(ns("loadings_bar"), height = "420px"),
               axes_note("the loading (weight) of each variable in the selected PC, from -1 to 1.",
                         "the original variables, sorted by loading.",
                         "Squared loadings in a PC add up to 1."))
        )
      ),

      # ==== 5. Reconstruction ==================================================
      nav_panel(
        "5 · Reconstruction", value = "recon",
        explain(withMathJax(
          p("PCA is also a way to ", strong("compress"), " data. Keep only each observation's first k",
            "scores, then run the recipe backwards to estimate its original values."),
          math_words(
            "$$\\hat{X} = T_k\\, V_k^\\top \\quad\\text{(then undo the scaling and add the averages back)}$$",
            "rebuild every observation from just its first k scores and the first k directions.",
            "With all the PCs you get the original data back exactly. With fewer, each point is pulled",
            "onto the best-fitting line (k = 1), plane (k = 2), and so on. The arrows in the chart show how",
            "far each point moves: these are the residuals from tab 1, now in more dimensions.",
            symbols = list(
              "\\(\\hat{X}\\)" = "the rebuilt data (the 'hat' means 'estimate')",
              "\\(T_k\\)" = "the first k columns of scores",
              "\\(V_k^\\top\\)" = "the first k directions, flipped so rows become columns"
            )
          )
        )),
        layout_columns(
          col_widths = c(4, 8),
          card(
            card_header("Controls"),
            sliderInput(ns("k"), "Components kept (k)", 1, 4, 1, step = 1),
            selectInput(ns("rx"), "Show variable (x)", NULL),
            selectInput(ns("ry"), "against variable (y)", NULL),
            uiOutput(ns("recon_stats"))
          ),
          card(full_screen = TRUE,
               card_header("Original (hollow) vs. reconstructed (filled)"),
               plotOutput(ns("recon_plot"), height = "440px"),
               axes_note("the variable chosen as x, in its original units.",
                         "the variable chosen as y, in its original units.",
                         "Each arrow goes from an observation's true values to its reconstruction from k PCs."))
        ),
        card(card_header("Reconstruction error by k"),
             plotOutput(ns("recon_curve"), height = "260px"),
             axes_note("k, the number of principal components kept.",
                       "% of the total variance that is lost (not captured by the first k PCs)."))
      ),

      # ==== 6. 3D ==============================================================
      nav_panel(
        "6 · Into 3D", value = "3d",
        explain(
          p("Here is a 3-D cloud that is mostly flat - like a tilted pancake.",
            "PC1 and PC2 span the plane of the pancake; PC3 points through its thin side.",
            "Dropping PC3 flattens the pancake onto the plane, and the right-hand plot shows",
            "that flat view. Drag the 3-D plot to rotate it, and make the cloud thicker to see",
            "how much information the flat view starts to lose.")
        ),
        layout_columns(
          col_widths = c(7, 5),
          card(full_screen = TRUE, card_header("3-D data with principal axes"),
               plotlyOutput(ns("plot3d"), height = "520px")),
          div(
            uiOutput(ns("stats3d")),
            card(card_header("The same data on PC1 vs PC2"),
                 plotOutput(ns("flat3d"), height = "360px"),
                 axes_note("each point's score on PC1 (position along the long side of the pancake).",
                           "each point's score on PC2 (position along the short side)."))
          )
        )
      ),

      # ==== 7. Quiz ============================================================
      nav_panel(
        "Quiz", value = "quiz",
        quiz_panel(ns)
      )
    )
  )
}

pca_server <- function(id) {
  moduleServer(id, function(input, output, session) {

    # ======================================================================
    # 1. Big idea
    # ======================================================================
    seed2 <- reactiveVal(1)
    observeEvent(input$resample2, seed2(sample.int(1e6, 1)))

    data2 <- reactive({
      X <- gen_2d(input$n2, input$rho, input$sdx, input$sdy, seed2())
      scale(X, scale = FALSE)
    })
    fit2 <- reactive(pca_fit(data2(), scale = FALSE))
    pc1_angle <- reactive({
      v <- fit2()$V[, 1]
      unname(atan2(v[2], v[1]) * 180 / pi) %% 180
    })
    revealed <- reactiveVal(FALSE)

    # Rotation is driven by rotate.js in the browser; the server just follows
    # the slider.
    angle <- reactive(input$angle %||% 20)
    observeEvent(input$snap, {
      revealed(TRUE)
      updateSliderInput(session, "angle", value = round(pc1_angle()))
    })
    observeEvent(list(input$n2, input$rho, input$sdx, input$sdy, input$resample2),
                 revealed(FALSE), ignoreInit = TRUE)

    output$idea_scatter <- renderPlot({
      X <- data2()
      th <- angle() * pi / 180
      u <- c(cos(th), sin(th))
      s <- drop(X %*% u)
      P <- outer(s, u)
      df <- data.frame(x = X[, 1], y = X[, 2], px = P[, 1], py = P[, 2])
      L <- max(abs(X)) * 1.15

      g <- ggplot(df) +
        geom_hline(yintercept = 0, colour = sv_col$grid) +
        geom_vline(xintercept = 0, colour = sv_col$grid)
      if (revealed()) {
        v <- fit2()$V
        g <- g +
          annotate("segment", x = -L * v[1, 1], y = -L * v[2, 1], xend = L * v[1, 1], yend = L * v[2, 1],
                   colour = sv_col$good, linewidth = 1, linetype = "22") +
          annotate("segment", x = -L * v[1, 2], y = -L * v[2, 2], xend = L * v[1, 2], yend = L * v[2, 2],
                   colour = sv_col$good, linewidth = 0.6, linetype = "22") +
          annotate("label", x = 0.85 * L * v[1, 1], y = 0.85 * L * v[2, 1], label = "PC1",
                   colour = sv_col$good, fill = sv_col$surface, label.size = 0, fontface = "bold") +
          annotate("label", x = 0.6 * L * v[1, 2], y = 0.6 * L * v[2, 2], label = "PC2",
                   colour = sv_col$good, fill = sv_col$surface, label.size = 0)
      }
      if (isTRUE(input$show_resid)) {
        g <- g + geom_segment(aes(x, y, xend = px, yend = py),
                              colour = sv_col$muted, linewidth = 0.35, alpha = 0.7)
      }
      g <- g + annotate("segment", x = -L * u[1], y = -L * u[2], xend = L * u[1], yend = L * u[2],
                        colour = sv_col$accent, linewidth = 1.1)
      if (isTRUE(input$show_proj)) {
        g <- g + geom_point(aes(px, py), colour = sv_col$accent, size = 1.8, alpha = 0.8)
      }
      # arc showing θ, measured counter-clockwise from the positive x1 axis
      r <- 0.5 * L
      arc_t <- seq(0, th, length.out = 80)
      g + geom_point(aes(x, y), colour = sv_col$primary, fill = sv_col$primary,
                     size = 2.4, alpha = 0.75, shape = 21, stroke = 0.5) +
        annotate("segment", x = 0, y = 0, xend = L, yend = 0,
                 colour = sv_col$ink2, linewidth = 0.6, linetype = "solid") +
        annotate("path", x = r * cos(arc_t), y = r * sin(arc_t),
                 colour = sv_col$accent, linewidth = 1.4,
                 arrow = arrow(length = unit(if (th > 0.1) 0.18 else 0, "cm"), type = "closed")) +
        annotate("label", x = r * cos(th / 2), y = r * sin(th / 2), label = "θ",
                 colour = sv_col$accent, fill = sv_col$surface, label.size = 0,
                 fontface = "bold", size = 5) +
        coord_equal(xlim = c(-L, L), ylim = c(-L, L), expand = FALSE) +
        labs(title = sprintf("θ = %d°", round(angle())),
             x = "x1  (variable 1, centered)", y = "x2  (variable 2, centered)") +
        theme_statviz()
    })

    output$idea_curve <- renderPlot({
      S <- fit2()$S
      m <- nrow(data2()) - 1   # SS = (n - 1) x variance
      tot <- m * sum(diag(S))
      ang <- 0:180
      kept <- m * vapply(ang, function(a) proj_variance(S, a), 0)
      df <- rbind(data.frame(angle = ang, value = kept, what = "Projected SS (kept)"),
                  data.frame(angle = ang, value = tot - kept, what = "Residual SS (lost)"))
      cur <- m * proj_variance(S, angle())
      g <- ggplot(df, aes(angle, value, colour = what)) +
        geom_hline(yintercept = tot, colour = sv_col$ink2, linetype = "dotted") +
        annotate("text", x = 2, y = tot, label = "total SS", vjust = -0.5, hjust = 0,
                 colour = sv_col$ink2, size = 3.3) +
        geom_line(linewidth = 1) +
        geom_vline(xintercept = angle(), colour = sv_col$accent, linewidth = 0.8) +
        annotate("point", x = angle(), y = cur, colour = sv_col$accent, size = 3.5) +
        annotate("point", x = angle(), y = tot - cur, colour = sv_col$muted, size = 3) +
        scale_colour_manual(values = c("Projected SS (kept)" = sv_col$primary,
                                       "Residual SS (lost)" = sv_col$muted), name = NULL) +
        scale_x_continuous(breaks = seq(0, 180, 45), limits = c(0, 180),
                           labels = function(x) paste0(x, "°")) +
        scale_y_continuous(limits = c(0, tot * 1.12)) +
        labs(x = "Angle θ of the line", y = "Sum of squares") +
        theme_statviz()
      if (revealed()) {
        g <- g + geom_vline(xintercept = pc1_angle(), colour = sv_col$good, linetype = "22") +
          annotate("label", x = pc1_angle(), y = tot * 0.06, label = "PC1",
                   colour = sv_col$good, fill = sv_col$surface, label.size = 0, fontface = "bold")
      }
      g
    })

    output$idea_stats <- renderUI({
      S <- fit2()$S
      m <- nrow(data2()) - 1
      tot <- m * sum(diag(S))
      kept <- m * proj_variance(S, angle())
      best <- m * fit2()$values[1]
      layout_column_wrap(
        width = 1 / 2, fill = FALSE,
        value_box(title = "Your line keeps",
                  value = sprintf("%.0f%%", 100 * kept / tot),
                  p(sprintf("SS %.0f of %.0f", kept, tot), br(),
                    sprintf("residual SS %.0f", tot - kept)),
                  theme = value_box_theme(bg = "#eaf2fc", fg = sv_col$ink)),
        value_box(title = if (revealed()) "PC1 keeps" else "Best possible",
                  value = if (revealed()) sprintf("%.0f%%", 100 * best / tot) else "?",
                  p(if (revealed()) tagList(sprintf("at %.0f°", pc1_angle()), br(),
                                            sprintf("residual SS %.0f (smallest possible)", tot - best))
                    else "press 'Reveal PC1'"),
                  theme = value_box_theme(bg = "#e6f6ef", fg = sv_col$ink))
      )
    })

    # ======================================================================
    # Real data (tabs 2-5)
    # ======================================================================
    raw <- reactive({
      if (input$dataset == "upload") {
        req(input$file)
        df <- tryCatch(read.csv(input$file$datapath, check.names = FALSE),
                       error = function(e) NULL)
        validate(need(!is.null(df), "Couldn't read that file as CSV."))
        list(name = input$file$name, data = df, group = NULL,
             blurb = sprintf("%d rows, %d columns.", nrow(df), ncol(df)))
      } else {
        pca_datasets[[input$dataset]]
      }
    })

    output$blurb <- renderUI(p(class = "small text-muted", raw()$blurb))

    observeEvent(raw(), {
      df <- raw()$data
      num <- names(df)[vapply(df, is.numeric, TRUE)]
      cat_cols <- names(df)[vapply(df, function(x) !is.numeric(x) && length(unique(x)) <= 8, TRUE)]
      updateSelectizeInput(session, "vars", choices = num, selected = num)
      updateSelectInput(session, "group", choices = c("(none)" = "", cat_cols),
                        selected = raw()$group %||% "")
    })

    clean <- reactive({
      df <- raw()$data
      vars <- intersect(input$vars, names(df))
      validate(need(length(vars) >= 2, "Pick at least two numeric variables."))
      keep <- complete.cases(df[, vars, drop = FALSE])
      df <- df[keep, , drop = FALSE]
      X <- as.matrix(df[, vars, drop = FALSE])
      sds <- apply(X, 2, sd)
      validate(need(all(sds > 0), "Drop variables that are constant - they have no variance to explain."))
      validate(need(nrow(X) > length(vars), "Need more rows than variables."))
      grp <- if (nzchar(input$group) && input$group %in% names(df)) factor(df[[input$group]]) else NULL
      labels <- if (input$dataset %in% c("USArrests", "mtcars", "swiss")) rownames(df) else
        paste("Row", which(keep))
      list(X = X, group = grp, labels = labels)
    })

    fit <- reactive(pca_fit(clean()$X, scale = input$scale))

    # keep PC pickers & reconstruction controls in sync with the data
    observeEvent(fit(), {
      pcs <- colnames(fit()$V)
      vars <- rownames(fit()$V)
      sel <- function(cur, default) if (isTRUE(cur %in% pcs)) cur else default
      updateSelectInput(session, "pcx", choices = pcs, selected = sel(input$pcx, pcs[1]))
      updateSelectInput(session, "pcy", choices = pcs, selected = sel(input$pcy, pcs[2]))
      updateSelectInput(session, "load_pc", choices = pcs, selected = sel(input$load_pc, pcs[1]))
      p <- length(pcs)
      updateSliderInput(session, "k", max = p, value = min(input$k, p))
      vsel <- function(cur, default) if (isTRUE(cur %in% vars)) cur else default
      updateSelectInput(session, "rx", choices = vars, selected = vsel(input$rx, vars[1]))
      updateSelectInput(session, "ry", choices = vars, selected = vsel(input$ry, vars[2]))
    })

    fmt_tbl <- function(M, digits = 3) {
      df <- as.data.frame(M)
      df[] <- lapply(df, function(x) if (is.numeric(x)) formatC(x, format = "f", digits = digits) else x)
      df
    }

    # ---- Step 1
    output$step1_summary <- renderTable({
      f <- fit()
      data.frame(Variable = names(f$means), Mean = f$means, SD = f$sds,
                 `Mean after` = colMeans(f$Z), `SD after` = apply(f$Z, 2, sd),
                 check.names = FALSE)
    }, digits = 2, striped = TRUE, spacing = "s")

    output$step1_head <- renderTable(head(fit()$Z, 6), digits = 2, spacing = "s")

    # ---- Step 2
    output$step2_cov <- renderPlot({
      mat_heatmap(fit()$S, if (input$scale) "correlation" else "covariance")
    })

    # ---- Step 3
    output$step3_eigs <- renderTable({
      f <- fit()
      data.frame(PC = names(f$values), `Eigenvalue (λ)` = f$values,
                 `Var. of scores` = apply(f$scores, 2, var),
                 `% of variance` = 100 * f$prop, check.names = FALSE)
    }, digits = 3, striped = TRUE, spacing = "s")

    output$step3_vecs <- renderPlot(mat_heatmap(fit()$V, "loading"))

    output$step3_check <- renderUI({
      f <- fit()
      pr <- prcomp(f$X, center = TRUE, scale. = f$scale)
      same <- isTRUE(all.equal(abs(unname(pr$rotation)), abs(unname(f$V)), tolerance = 1e-6))
      div(class = if (same) "alert alert-success py-2 small" else "alert alert-warning py-2 small",
          icon(if (same) "check" else "triangle-exclamation"),
          if (same) "Sanity check: R's built-in prcomp() gives the same eigenvectors (up to sign)."
          else "prcomp() disagrees here - probably because of repeated eigenvalues, where the axes aren't unique.")
    })

    # ---- Step 3: covariance -> eigenvector bridge (uses the tab-1 cloud)
    u_vec <- reactive({
      th <- (input$u_angle %||% 110) * pi / 180
      c(cos(th), sin(th))
    })

    observeEvent(input$power_step, {
      Su <- unname(drop(fit2()$S %*% u_vec()))
      updateSliderInput(session, "u_angle", value = round((atan2(Su[2], Su[1]) * 180 / pi) %% 180))
    })

    output$bridge_plot <- renderPlot({
      X <- data2(); S <- fit2()$S; V <- fit2()$V
      u <- u_vec(); Su <- drop(S %*% u)
      t <- seq(0, 2 * pi, length.out = 200)
      circ <- cbind(cos(t), sin(t))
      ell <- circ %*% S            # S is symmetric
      L <- max(c(abs(ell), abs(Su), 1.2)) * 1.25
      ev <- data.frame(x = -L * V[1, ], y = -L * V[2, ], xend = L * V[1, ], yend = L * V[2, ],
                       lab = paste0("v", 1:2))
      ggplot() +
        geom_point(data = data.frame(x = X[, 1], y = X[, 2]), aes(x, y),
                   colour = sv_col$muted, alpha = 0.3, size = 1.6) +
        geom_segment(data = ev, aes(x, y, xend = xend, yend = yend),
                     colour = sv_col$good, linetype = "22", linewidth = 0.7) +
        geom_text(data = ev, aes(0.92 * xend, 0.92 * yend, label = lab), colour = sv_col$good,
                  fontface = "bold", vjust = -0.6) +
        geom_path(data = data.frame(x = circ[, 1], y = circ[, 2]), aes(x, y),
                  colour = sv_col$ink2, linewidth = 0.5) +
        geom_path(data = data.frame(x = ell[, 1], y = ell[, 2]), aes(x, y),
                  colour = sv_col$primary, linewidth = 0.8, alpha = 0.6) +
        annotate("segment", x = 0, y = 0, xend = Su[1], yend = Su[2], colour = sv_col$primary,
                 linewidth = 1.6, arrow = arrow(length = unit(0.3, "cm"), type = "closed")) +
        annotate("segment", x = 0, y = 0, xend = u[1], yend = u[2], colour = sv_col$accent,
                 linewidth = 1.6, arrow = arrow(length = unit(0.3, "cm"), type = "closed")) +
        annotate("label", x = 1.15 * u[1], y = 1.15 * u[2], label = "u", colour = sv_col$accent,
                 fill = sv_col$surface, label.size = 0, fontface = "bold", size = 5) +
        annotate("label", x = Su[1] * 1.08 + 0.12 * sign(Su[1]), y = Su[2] * 1.08, label = "Su",
                 colour = sv_col$primary, fill = sv_col$surface, label.size = 0, fontface = "bold", size = 5) +
        coord_equal(xlim = c(-L, L), ylim = c(-L, L), expand = FALSE) +
        labs(x = "x1 (centered)", y = "x2 (centered)") +
        theme_statviz()
    })

    output$bridge_math <- renderUI({
      S <- fit2()$S; u <- u_vec(); Su <- drop(S %*% u)
      gap <- acos(min(1, abs(sum(u * Su)) / sqrt(sum(Su^2)))) * 180 / pi
      f2 <- function(x) formatC(x, format = "f", digits = 2)
      mat <- sprintf("\\begin{pmatrix} %s & %s \\\\ %s & %s \\end{pmatrix}",
                     f2(S[1, 1]), f2(S[1, 2]), f2(S[2, 1]), f2(S[2, 2]))
      vec <- function(v) sprintf("\\begin{pmatrix} %s \\\\ %s \\end{pmatrix}", f2(v[1]), f2(v[2]))
      eig <- gap < 1
      withMathJax(tagList(
        p(class = "mb-1", sprintf("$$\\underset{S}{%s}\\,\\underset{u}{%s} = \\underset{Su}{%s}$$",
                                  mat, vec(u), vec(Su))),
        p(class = "small text-muted",
          "Each entry of Su: go along a row of S, multiply by the matching entry of u, and add."),
        p(class = "mb-1", "Spread (variance) of the data along u: \\(u^\\top S u =\\)", strong(f2(sum(u * Su)))),
        p(class = "mb-1", "Angle between u and Su:", strong(sprintf("%.1f°", gap)),
          span(class = "text-muted small", "(0° means S didn't turn u)")),
        if (eig) {
          lam <- sum(u * Su)
          div(class = "alert alert-success py-2 small mt-2", icon("check"),
              sprintf(" u is an eigenvector! S only stretched it, by λ = %s.", f2(lam)),
              if (lam >= mean(fit2()$values))
                " That stretch factor is the spread along u, the largest possible: u is PC1."
              else
                " That stretch factor is the spread along u, the smallest possible: u is PC2.")
        } else {
          div(class = "alert alert-secondary py-2 small mt-2",
              "S turned the arrow, so u is not an eigenvector. Keep turning u, or click the button.")
        }
      ))
    })

    output$step3_verify <- renderTable({
      f <- fit()
      k <- min(2, ncol(f$V))
      out <- data.frame(Variable = rownames(f$V))
      for (j in seq_len(k)) {
        out[[sprintf("S × v%d", j)]] <- drop(f$S %*% f$V[, j])
        out[[sprintf("λ%d × v%d", j, j)]] <- f$values[j] * f$V[, j]
      }
      out
    }, digits = 3, striped = TRUE, spacing = "s")

    # ---- Step 5
    output$step5_example <- renderUI({
      f <- fit()
      w <- f$V[, 1]
      terms <- sprintf("%s %.2f × %s", ifelse(w < 0, "−", "+"), abs(w), names(w))
      terms[1] <- sub("^\\+ ", "", terms[1])
      first <- f$Z[1, ]
      div(class = "sv-math small",
          p(class = "mt-2 mb-1", strong("Example: "), "the recipe for PC1 on this dataset is"),
          p(class = "font-monospace mb-1",
            "score on PC1 = ", paste(terms, collapse = " ")),
          p(class = "text-muted mb-0",
            if (f$scale) "(using each variable's centered and scaled value)" else "(using each variable's centered value)",
            sprintf("For the first observation that gives %.2f.", sum(w * first))))
    })

    output$step4_scores <- renderPlotly({
      f <- fit(); cl <- clean()
      req(input$pcx %in% colnames(f$scores), input$pcy %in% colnames(f$scores))
      df <- data.frame(x = f$scores[, input$pcx], y = f$scores[, input$pcy], label = cl$labels)
      pct <- function(pc) sprintf("%s (%.1f%%)", pc, 100 * f$prop[pc])
      hover <- sprintf("<b>%s</b><br>%s: %.2f<br>%s: %.2f", df$label, input$pcx, df$x, input$pcy, df$y)
      if (!is.null(cl$group)) {
        df$group <- cl$group
        hover <- paste0(hover, "<br>", input$group, ": ", df$group)
        p <- plot_ly(df, x = ~x, y = ~y, color = ~group, colors = group_palette(levels(df$group)),
                     text = hover, hoverinfo = "text", type = "scatter", mode = "markers",
                     marker = list(size = 9, opacity = 0.8, line = list(color = sv_col$surface, width = 1.5)))
      } else {
        p <- plot_ly(df, x = ~x, y = ~y, text = hover, hoverinfo = "text", type = "scatter",
                     mode = "markers",
                     marker = list(size = 9, opacity = 0.8, color = sv_col$primary,
                                   line = list(color = sv_col$surface, width = 1.5)))
      }
      p |> layout(xaxis = list(title = pct(input$pcx), zerolinecolor = sv_col$grid, gridcolor = sv_col$grid),
                  yaxis = list(title = pct(input$pcy), zerolinecolor = sv_col$grid, gridcolor = sv_col$grid,
                               scaleanchor = "x"),
                  legend = list(orientation = "h", y = 1.08),
                  paper_bgcolor = sv_col$surface, plot_bgcolor = sv_col$surface) |>
        config(displaylogo = FALSE)
    })

    # ---- Scree
    output$scree_plot <- renderPlot({
      f <- fit()
      df <- data.frame(pc = factor(names(f$values), levels = names(f$values)),
                       prop = 100 * f$prop, cum = 100 * cumsum(f$prop))
      need_k <- which(df$cum >= input$thresh)[1]
      df$kept <- seq_len(nrow(df)) <= need_k
      ggplot(df, aes(pc)) +
        geom_col(aes(y = prop, fill = kept), width = 0.7) +
        geom_hline(yintercept = input$thresh, colour = sv_col$accent, linetype = "22") +
        annotate("text", x = 0.5, y = input$thresh, label = paste0(input$thresh, "% target"),
                 hjust = 0, vjust = -0.5, colour = sv_col$accent, size = 3.5) +
        geom_line(aes(y = cum, group = 1), colour = sv_col$ink2, linewidth = 0.8) +
        geom_point(aes(y = cum), colour = sv_col$ink2, size = 2.5) +
        geom_text(aes(y = prop, label = sprintf("%.0f%%", prop),
                      vjust = ifelse(prop > 12, 1.6, -0.5),
                      colour = I(ifelse(prop > 12 & kept, "#ffffff", sv_col$ink))),
                  size = 3.5, fontface = "bold") +
        annotate("text", x = nrow(df), y = 100, label = "cumulative", hjust = 1, vjust = -0.8,
                 colour = sv_col$ink2, size = 3.5) +
        scale_fill_manual(values = c(`TRUE` = sv_col$primary, `FALSE` = "#b7d3f6"),
                          labels = c(`TRUE` = "kept", `FALSE` = "dropped"), name = NULL) +
        scale_y_continuous(limits = c(0, 106), labels = function(x) paste0(x, "%")) +
        labs(x = NULL, y = "% of total variance") +
        theme_statviz()
    })

    output$scree_answer <- renderUI({
      f <- fit()
      cum <- cumsum(f$prop)
      k <- which(cum >= input$thresh / 100)[1]
      kaiser <- sum(f$values > 1)
      tagList(
        p("Keep ", strong(k), if (k == 1) " component" else " components",
          sprintf(" (%.1f%% of variance, down from %d variables).", 100 * cum[k], length(cum))),
        if (input$scale) p(class = "small text-muted",
                           "Kaiser rule says keep ", strong(kaiser), ".")
        else p(class = "small text-muted", "Kaiser rule only applies to scaled data.")
      )
    })

    output$scree_table <- renderTable({
      f <- fit()
      data.frame(PC = names(f$values), λ = f$values, `%` = 100 * f$prop,
                 `Cum. %` = 100 * cumsum(f$prop), check.names = FALSE)
    }, digits = 2, striped = TRUE, spacing = "s")

    # ---- Interpreting
    output$biplot <- renderPlot({
      f <- fit(); cl <- clean()
      req(input$pcx %in% colnames(f$scores), input$pcy %in% colnames(f$scores))
      sc <- data.frame(x = f$scores[, input$pcx], y = f$scores[, input$pcy])
      ld <- data.frame(var = rownames(f$V), x = f$V[, input$pcx], y = f$V[, input$pcy])
      # stretch arrows to fill the score cloud
      k <- 0.8 * max(abs(as.matrix(sc))) / max(sqrt(ld$x^2 + ld$y^2))
      ld$x <- ld$x * k; ld$y <- ld$y * k
      pct <- function(pc) sprintf("%s (%.1f%%)", pc, 100 * f$prop[pc])

      g <- ggplot() +
        geom_hline(yintercept = 0, colour = sv_col$grid) +
        geom_vline(xintercept = 0, colour = sv_col$grid)
      if (!is.null(cl$group)) {
        sc$group <- cl$group
        g <- g + geom_point(data = sc, aes(x, y, colour = group), size = 2.3, alpha = 0.7) +
          scale_colour_manual(values = group_palette(levels(sc$group)), name = input$group)
      } else {
        g <- g + geom_point(data = sc, aes(x, y), colour = sv_col$primary, size = 2.3, alpha = 0.6)
      }
      g + geom_segment(data = ld, aes(0, 0, xend = x, yend = y),
                       arrow = arrow(length = unit(0.25, "cm")),
                       colour = sv_col$ink, linewidth = 0.7) +
        ggrepel::geom_label_repel(data = ld, aes(x, y, label = var), colour = sv_col$ink,
                                  fill = scales::alpha(sv_col$surface, 0.85), label.size = 0,
                                  size = 3.8, min.segment.length = Inf, seed = 1) +
        coord_equal() +
        labs(x = pct(input$pcx), y = pct(input$pcy)) +
        theme_statviz()
    })

    output$loadings_bar <- renderPlot({
      f <- fit()
      req(input$load_pc %in% colnames(f$V))
      df <- data.frame(var = rownames(f$V), w = f$V[, input$load_pc])
      df$var <- factor(df$var, levels = df$var[order(df$w)])
      lim <- max(abs(df$w))
      ggplot(df, aes(w, var, fill = w)) +
        geom_col(width = 0.7) +
        geom_vline(xintercept = 0, colour = sv_col$ink2) +
        geom_text(aes(label = sprintf("%.2f", w), hjust = ifelse(w >= 0, -0.15, 1.15)),
                  colour = sv_col$ink, size = 3.6) +
        scale_fill_diverging(lim) +
        scale_x_continuous(limits = c(-lim, lim) * 1.3) +
        labs(x = sprintf("Loading on %s", input$load_pc), y = NULL,
             subtitle = sprintf("%s explains %.1f%% of variance", input$load_pc,
                                100 * f$prop[input$load_pc])) +
        theme_statviz() + theme(legend.position = "none")
    })

    # ---- Reconstruction
    recon <- reactive({
      f <- fit()
      k <- min(input$k, ncol(f$V))
      list(k = k, Xhat = pca_reconstruct(f, k))
    })

    output$recon_plot <- renderPlot({
      f <- fit(); r <- recon(); cl <- clean()
      req(input$rx %in% colnames(f$X), input$ry %in% colnames(f$X))
      df <- data.frame(x = f$X[, input$rx], y = f$X[, input$ry],
                       hx = r$Xhat[, input$rx], hy = r$Xhat[, input$ry])
      ggplot(df) +
        geom_segment(aes(x, y, xend = hx, yend = hy), colour = sv_col$muted, linewidth = 0.4,
                     arrow = arrow(length = unit(0.12, "cm"))) +
        geom_point(aes(x, y), shape = 21, colour = sv_col$primary, fill = NA, size = 2.6, stroke = 0.8) +
        geom_point(aes(hx, hy), colour = sv_col$accent, size = 2.2) +
        labs(x = input$rx, y = input$ry,
             subtitle = sprintf("Using %d of %d components", r$k, ncol(f$V))) +
        theme_statviz()
    })

    output$recon_stats <- renderUI({
      f <- fit(); r <- recon()
      kept <- 100 * sum(f$prop[seq_len(r$k)])
      p(class = "mt-2", "These ", strong(r$k), " PCs retain ", strong(sprintf("%.1f%%", kept)),
        " of the ", if (f$scale) "(scaled) ", "variance. The rest is the reconstruction error.")
    })

    output$recon_curve <- renderPlot({
      f <- fit()
      p <- ncol(f$V)
      lost <- 100 * (1 - cumsum(f$prop))
      df <- data.frame(k = seq_len(p), lost = pmax(lost, 0))
      ggplot(df, aes(k, lost)) +
        geom_line(colour = sv_col$primary, linewidth = 1) +
        geom_point(colour = sv_col$primary, size = 3) +
        annotate("point", x = recon()$k, y = df$lost[recon()$k], colour = sv_col$accent, size = 5) +
        scale_x_continuous(breaks = seq_len(p)) +
        scale_y_continuous(limits = c(0, NA), labels = function(x) paste0(x, "%")) +
        labs(x = "Components kept (k)", y = "Variance lost") +
        theme_statviz()
    })

    # ======================================================================
    # 6. 3D
    # ======================================================================
    seed3 <- reactiveVal(7)
    observeEvent(input$resample3, seed3(sample.int(1e6, 1)))
    data3 <- reactive(gen_3d(input$n3, input$thick, seed3()))
    fit3 <- reactive(pca_fit(data3(), scale = FALSE))

    output$plot3d <- renderPlotly({
      X <- data3(); f <- fit3()
      m <- f$means; V <- f$V; sd <- sqrt(f$values)
      p <- plot_ly() |>
        add_trace(x = X[, 1], y = X[, 2], z = X[, 3], type = "scatter3d", mode = "markers",
                  marker = list(size = 3, color = sv_col$primary, opacity = 0.75),
                  name = "data", hoverinfo = "skip")
      if (isTRUE(input$show_plane)) {
        a <- 2.5 * sd[1]; b <- 2.5 * sd[2]
        corners <- rbind(m + a * V[, 1] + b * V[, 2], m + a * V[, 1] - b * V[, 2],
                         m - a * V[, 1] - b * V[, 2], m - a * V[, 1] + b * V[, 2])
        p <- p |> add_trace(x = corners[, 1], y = corners[, 2], z = corners[, 3],
                            i = c(0, 0), j = c(1, 2), k = c(2, 3), type = "mesh3d",
                            color = I(sv_col$good), opacity = 0.25, name = "PC1-PC2 plane",
                            hoverinfo = "skip")
      }
      if (isTRUE(input$show_axes)) {
        cols <- c(sv_col$good, sv_col$good, sv_col$accent)
        for (j in 1:3) {
          end <- m + 2 * max(sd[j], 0.3) * V[, j]
          p <- p |> add_trace(x = c(m[1], end[1]), y = c(m[2], end[2]), z = c(m[3], end[3]),
                              type = "scatter3d", mode = "lines+text",
                              line = list(color = cols[j], width = 7),
                              text = c("", paste0("PC", j)), textfont = list(color = cols[j], size = 14),
                              name = paste0("PC", j), hoverinfo = "name")
        }
      }
      p |> layout(scene = list(aspectmode = "data",
                               xaxis = list(title = "x1"), yaxis = list(title = "x2"),
                               zaxis = list(title = "x3")),
                  showlegend = FALSE, paper_bgcolor = sv_col$surface,
                  margin = list(l = 0, r = 0, t = 0, b = 0)) |>
        config(displaylogo = FALSE)
    })

    output$stats3d <- renderUI({
      pr <- 100 * fit3()$prop
      layout_column_wrap(
        width = 1 / 3, fill = FALSE,
        value_box("PC1", sprintf("%.0f%%", pr[1]), theme = value_box_theme(bg = "#e6f6ef", fg = sv_col$ink)),
        value_box("PC2", sprintf("%.0f%%", pr[2]), theme = value_box_theme(bg = "#e6f6ef", fg = sv_col$ink)),
        value_box("PC3 (dropped)", sprintf("%.0f%%", pr[3]),
                  theme = value_box_theme(bg = "#fdeee7", fg = sv_col$ink))
      )
    })

    output$flat3d <- renderPlot({
      sc <- fit3()$scores
      ggplot(data.frame(x = sc[, 1], y = sc[, 2]), aes(x, y)) +
        geom_hline(yintercept = 0, colour = sv_col$grid) +
        geom_vline(xintercept = 0, colour = sv_col$grid) +
        geom_point(colour = sv_col$primary, alpha = 0.7, size = 2) +
        coord_equal() +
        labs(x = "PC1", y = "PC2") +
        theme_statviz()
    })

    # ======================================================================
    # 7. Quiz
    # ======================================================================
    quiz <- list(
      list(q = "PC1 is the direction that...",
           choices = c("passes through the most data points",
                       "maximizes the variance of the projected data",
                       "is parallel to the variable with the largest mean",
                       "minimizes the vertical (y) distance to the points, like regression"),
           answer = 2,
           why = "PC1 maximizes projected variance - equivalently it minimizes the perpendicular (not vertical) distances."),
      list(q = "Why do we usually scale variables before PCA?",
           choices = c("It makes the eigenvalues add up to 1",
                       "Otherwise variables measured in large units dominate the components",
                       "PCA only works on normally distributed data",
                       "It removes outliers"),
           answer = 2,
           why = "Unscaled PCA uses the covariance matrix, so a variable with variance 1000 swamps one with variance 1. Try USArrests with and without scaling."),
      list(q = "An eigenvalue in PCA represents...",
           choices = c("the angle of the component",
                       "the correlation between two variables",
                       "the variance of the data along that component",
                       "the number of observations in a cluster"),
           answer = 3,
           why = "λ_k is the variance of the scores on PC k; λ_k / Σλ is its share of the total."),
      list(q = "The scores on PC1 and PC2 are...",
           choices = c("perfectly correlated", "uncorrelated",
                       "negatively correlated", "it depends on the dataset"),
           answer = 2,
           why = "The eigenvectors of a symmetric matrix are orthogonal, so the projected scores have zero correlation."),
      list(q = "Your loadings for PC1 are all positive and roughly equal. PC1 is best described as...",
           choices = c("a contrast between variables", "an overall 'size' or average of the variables",
                       "noise", "the variable with the largest loading"),
           answer = 2,
           why = "Equal positive weights make PC1 essentially a (scaled) average - common in body-measurement data like iris.")
    )

    quiz_server(input, output, session, quiz)
  })
}
