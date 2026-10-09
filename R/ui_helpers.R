# Shared UI building blocks for every StatViz module ---------------------------

explain <- function(...) {
  div(class = "sv-explain", ...)
}

# A formula followed by a plain-language reading of it, and optionally a key
# to its symbols: math_words("$$...$$", "what it says", symbols = list("\\(x\\)" = "..."))
math_words <- function(formula, ..., symbols = NULL) {
  div(
    class = "sv-math",
    div(class = "sv-formula", formula),
    div(class = "sv-words", strong("In words: "), ...),
    if (length(symbols)) div(
      class = "sv-symbols",
      strong("Symbols: "),
      lapply(names(symbols), function(sym) span(class = "sv-sym", sym, " ", symbols[[sym]]))
    )
  )
}

# Extra math for the curious, collapsed by default.
deeper <- function(summary, ...) {
  tags$details(class = "sv-deeper", tags$summary(icon("calculator"), summary), div(...))
}

# "How to read this plot" footer for a card.
axes_note <- function(x, y, ...) {
  card_footer(
    class = "small sv-axes",
    div(strong("Horizontal axis: "), x),
    if (nzchar(y)) div(strong("Vertical axis: "), y),
    if (length(list(...))) div(class = "text-muted mt-1", ...)
  )
}

# Multiple-choice quiz ---------------------------------------------------------
# questions: list of list(q = "...", choices = c(...), answer = <index>, why = "...")

quiz_panel <- function(ns) {
  card(
    card_header("Check your understanding"),
    uiOutput(ns("quiz_questions")),
    actionButton(ns("quiz_check"), "Check my answers", class = "btn-primary", width = "220px"),
    uiOutput(ns("quiz_score"))
  )
}

quiz_server <- function(input, output, session, questions) {
  output$quiz_questions <- renderUI({
    ns <- session$ns
    tagList(lapply(seq_along(questions), function(i) {
      q <- questions[[i]]
      div(class = "mb-3",
          radioButtons(ns(paste0("q", i)), paste0(i, ". ", q$q),
                       choiceNames = q$choices, choiceValues = seq_along(q$choices),
                       selected = character(0), width = "100%"),
          uiOutput(ns(paste0("fb", i))))
    }))
  })

  observeEvent(input$quiz_check, {
    correct <- 0
    for (i in seq_along(questions)) local({
      ii <- i
      ans <- input[[paste0("q", ii)]]
      right <- !is.null(ans) && as.integer(ans) == questions[[ii]]$answer
      if (right) correct <<- correct + 1
      output[[paste0("fb", ii)]] <- renderUI({
        if (is.null(ans)) {
          div(class = "small text-muted", "No answer yet.")
        } else {
          div(class = if (right) "alert alert-success py-1 small" else "alert alert-danger py-1 small",
              icon(if (right) "check" else "xmark"),
              strong(if (right) " Correct. " else " Not quite. "), questions[[ii]]$why)
        }
      })
    })
    output$quiz_score <- renderUI(
      h5(class = "mt-3", sprintf("Score: %d / %d", correct, length(questions)))
    )
  })
}

# Sidebar section that shows only on some tabs of a module's navset.
show_on <- function(ns, tabs, ...) {
  conditionalPanel(sprintf("[%s].includes(input.tab)", paste0("'", tabs, "'", collapse = ",")),
                   ns = ns, ...)
}

# Value box in the app's soft palette.
sv_box <- function(title, value, ..., tone = c("blue", "green", "orange", "gray")) {
  bg <- c(blue = "#eaf2fc", green = "#e6f6ef", orange = "#fdeee7", gray = "#f0efec")[match.arg(tone)]
  value_box(title = title, value = value, ..., theme = value_box_theme(bg = bg, fg = sv_col$ink))
}
