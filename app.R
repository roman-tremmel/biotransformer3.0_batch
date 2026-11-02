# =============================================================================
# BioTransformer 3.0 - batch metabolism prediction (Shiny front-end)
#
# Paste a list of SMILES (or upload a CSV) and predict small-molecule
# metabolism for the whole set at once. Each input keeps its own row even when
# it cannot be parsed or predicted, and the batch always runs to completion.
#
# The heavy lifting is done by the official BioTransformer 3.0 JAR via
# R/biotransformer.R. Run setup.sh once to fetch and prepare the JAR.
# =============================================================================

library(shiny)
library(DT)

app_dir <- normalizePath(".", mustWork = FALSE)
source(file.path("R", "biotransformer.R"))

# ---- input parsing ---------------------------------------------------------

# Accept free-text (one entry per line) in any of:
#   SMILES
#   name,SMILES   /   name;SMILES   /   name<TAB>SMILES   /   name<space>SMILES
# A SMILES token is detected as the field that looks most like a SMILES.
.looks_like_smiles <- function(x) grepl("[=#\\[\\]()@+]|[A-Za-z].*[0-9]|^[A-Za-z]{2,}$", x) &&
  grepl("^[^ ]+$", x)

parse_smiles_text <- function(txt) {
  lines <- unlist(strsplit(txt, "\r?\n"))
  lines <- trimws(lines)
  lines <- lines[nzchar(lines) & !grepl("^#", lines)]
  if (!length(lines)) return(data.frame(name = character(0), smiles = character(0)))
  # drop an obvious header line
  if (grepl("^(name|id|compound)?[,;\t ]*smiles$", lines[1], ignore.case = TRUE) ||
      tolower(lines[1]) == "smiles")
    lines <- lines[-1]
  name <- character(length(lines)); smi <- character(length(lines))
  for (i in seq_along(lines)) {
    parts <- unlist(strsplit(lines[i], "[,;\t]| {2,}|\\s+"))
    parts <- parts[nzchar(parts)]
    if (length(parts) >= 2) {
      # pick the token with the most SMILES-ish characters as the structure
      score <- nchar(gsub("[^=#\\[\\]()@+0-9]", "", parts))
      j <- which.max(score); smi[i] <- parts[j]; name[i] <- paste(parts[-j], collapse = " ")
    } else {
      smi[i] <- parts[1]; name[i] <- ""
    }
  }
  data.frame(name = name, smiles = smi, stringsAsFactors = FALSE)
}

parse_csv_file <- function(path) {
  df <- tryCatch(utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
                 error = function(e) NULL)
  if (is.null(df) || !ncol(df)) return(NULL)
  nm <- tolower(trimws(names(df)))
  scol <- which(nm %in% c("smiles", "smile", "structure", "canonical_smiles"))
  ncol_ <- which(nm %in% c("name", "id", "compound", "compound_name", "title"))
  if (!length(scol)) {
    # file may be headerless like the vendor example: id,name,smiles,...
    df2 <- tryCatch(utils::read.csv(path, stringsAsFactors = FALSE, header = FALSE),
                    error = function(e) NULL)
    if (!is.null(df2) && ncol(df2) >= 3) {
      return(data.frame(name = as.character(df2[[2]]),
                        smiles = as.character(df2[[3]]), stringsAsFactors = FALSE))
    }
    # otherwise assume first column is SMILES
    return(data.frame(name = "", smiles = as.character(df[[1]]), stringsAsFactors = FALSE))
  }
  data.frame(name = if (length(ncol_)) as.character(df[[ncol_[1]]]) else "",
             smiles = as.character(df[[scol[1]]]), stringsAsFactors = FALSE)
}

# ---- UI --------------------------------------------------------------------

badge <- function(txt, cls) sprintf('<span class="bt-badge %s">%s</span>', cls, txt)
status_class <- c("predicted" = "bt-ok", "no metabolites" = "bt-none",
                  "parse error" = "bt-parse", "error" = "bt-err",
                  "timeout" = "bt-err")

ui <- fluidPage(
  tags$head(
    tags$title("BioTransformer 3.0 - batch"),
    tags$style(HTML("
      body { background:#f7f8fa; }
      .bt-title { font-weight:700; }
      .bt-sub { color:#5b6370; margin-top:-6px; }
      .bt-badge { display:inline-block; padding:2px 8px; border-radius:10px;
                  font-size:12px; font-weight:600; margin-right:4px; }
      .bt-ok    { background:#e6f4ea; color:#137333; }
      .bt-none  { background:#eef1f4; color:#5b6370; }
      .bt-parse { background:#fff4e5; color:#9a6400; }
      .bt-err   { background:#fdeaea; color:#b3261e; }
      .bt-legend { font-size:12.5px; color:#5b6370; margin:6px 0 2px; }
      .bt-struct { width:120px; height:84px; }
      .bt-tiles .shiny-html-output { display:inline-block; }
      .bt-tile { display:inline-block; background:#fff; border:1px solid #e3e6ea;
                 border-radius:8px; padding:6px 12px; margin:4px 6px 0 0; }
      .bt-tile b { font-size:18px; }
      .well { background:#fff; }
    ")),
    tags$script(src = "smiles-drawer.min.js"),
    tags$script(HTML("
      // redraw all metabolite structure cells currently in the DOM
      window.btDrawStructures = function(){
        try {
          if (typeof SmilesDrawer === 'undefined') return;
          SmilesDrawer.apply({width:120, height:84, compactDrawing:false, padding:2},
                             {}, 'data-smiles', 'light');
        } catch(e) {}
      };
      // single-molecule live preview
      Shiny.addCustomMessageHandler('bt_preview', function(msg){
        try {
          var el = document.getElementById('bt_preview_svg');
          if (!el) return;
          if (!msg.smiles) { el.innerHTML=''; return; }
          if (typeof SmiDrawer === 'undefined') return;
          var d = new SmiDrawer({width:260, height:170}, {});
          d.draw(msg.smiles, '#bt_preview_svg', 'light', null, function(){
            el.innerHTML='<text x=\"8\" y=\"20\" fill=\"#b3261e\" font-size=\"12\">invalid SMILES</text>';
          });
        } catch(e) {}
      });
    "))
  ),

  titlePanel(div(
    div(class = "bt-title", "BioTransformer 3.0 — batch metabolism prediction"),
    div(class = "bt-sub",
        "Predict in-silico metabolites for many compounds at once. Paste SMILES or upload a CSV.")
  )),

  sidebarLayout(
    sidebarPanel(
      width = 4,
      uiOutput("install_warning"),
      tags$label("SMILES input", style = "font-weight:600;"),
      helpText("One compound per line. Optionally prefix a name: ",
               tags$code("name, SMILES"), "."),
      textAreaInput("smiles_text", NULL, rows = 7,
                    placeholder = "caffeine, CN1C=NC2=C1C(=O)N(C(=O)N2C)C\nthymol, CC(C)C1=CC=C(C)C=C1O"),
      fileInput("smiles_file", "… or upload a CSV (a 'SMILES' column, optional 'name')",
                accept = c(".csv", ".txt")),
      actionLink("load_example", "load example set"),
      tags$hr(),
      selectInput("btType", "Biotransformer", choices = BT_TYPES, selected = "allHuman"),
      conditionalPanel(
        "input.btType == 'allHuman' || input.btType == 'cyp450'",
        selectInput("cyp_mode", "CYP450 Mode", choices = BT_CYP_MODES, selected = "3"),
        helpText(HTML("<b>Combined</b>: CypReact + rules + CyProduct &middot; ",
                      "<b>Rule Based</b>: CypReact + BioTransformer rules &middot; ",
                      "<b>CyProduct</b>: CyProduct only"))),
      sliderInput("steps", "Prediction steps", min = 1, max = 3, value = 1, step = 1),
      fluidRow(
        column(6, numericInput("workers", "Parallel jobs", value = 2, min = 1, max = 16, step = 1)),
        column(6, numericInput("timeout", "Per-molecule limit (s)", value = 300, min = 30, step = 30))
      ),
      checkboxInput("offline", "Offline mode (skip online PubChem/synonym annotation — faster)", value = FALSE),
      checkboxInput("show_struct", "Show metabolite structures", value = TRUE),
      tags$hr(),
      actionButton("run", "Run prediction", class = "btn-primary", width = "100%"),
      tags$div(style = "margin-top:10px;",
               tags$label("Structure preview"),
               tags$div(style = "background:#fff;border:1px solid #e3e6ea;border-radius:8px;padding:4px;",
                        HTML('<svg id="bt_preview_svg" width="260" height="170"></svg>')))
    ),

    mainPanel(
      width = 8,
      uiOutput("status_tiles"),
      conditionalPanel("output.has_results == true",
        div(class = "bt-legend",
            HTML(paste("Status:",
                       badge("predicted", "bt-ok"), "metabolites found &middot;",
                       badge("no metabolites", "bt-none"), "parsed, none predicted &middot;",
                       badge("parse error", "bt-parse"), "SMILES unreadable &middot;",
                       badge("error / timeout", "bt-err"), "skipped, batch continued")))
      ),
      tabsetPanel(
        id = "tabs",
        tabPanel("Metabolites",
          br(),
          fluidRow(
            column(6, downloadButton("dl_results", "Download metabolites (CSV)")),
            column(6, div(style = "text-align:right;", textOutput("result_count", inline = TRUE)))
          ),
          br(), DT::dataTableOutput("results_tbl")
        ),
        tabPanel("Per-compound summary",
          br(),
          downloadButton("dl_summary", "Download summary (CSV)"),
          br(), br(), DT::dataTableOutput("summary_tbl")
        )
      )
    )
  )
)

# ---- server ----------------------------------------------------------------

server <- function(input, output, session) {

  problems <- bt_check_install(app_dir)
  output$install_warning <- renderUI({
    if (!length(problems)) return(NULL)
    div(style = "background:#fdeaea;border:1px solid #f3b4ae;border-radius:8px;padding:8px 10px;margin-bottom:10px;color:#8a1c12;",
        tags$b("Setup incomplete:"), tags$ul(lapply(problems, tags$li)),
        "Run ", tags$code("bash setup.sh"), " first.")
  })

  rv <- reactiveValues(results = NULL, summary = NULL)

  # live single-molecule preview: use the first input line
  observeEvent(input$smiles_text, {
    first <- parse_smiles_text(input$smiles_text)
    smi <- if (nrow(first)) first$smiles[1] else ""
    session$sendCustomMessage("bt_preview", list(smiles = smi))
  }, ignoreInit = TRUE)

  observeEvent(input$load_example, {
    ex <- "caffeine, CN1C=NC2=C1C(=O)N(C(=O)N2C)C\nthymol, CC(C)C1=CC=C(C)C=C1O\natrazine, CCNC1=NC(=NC(=N1)Cl)NC(C)C\nacetaminophen, CC(=O)NC1=CC=C(O)C=C1"
    updateTextAreaInput(session, "smiles_text", value = ex)
  })

  current_input <- reactive({
    df <- NULL
    if (!is.null(input$smiles_file)) df <- parse_csv_file(input$smiles_file$datapath)
    if (is.null(df) || !nrow(df)) df <- parse_smiles_text(input$smiles_text %||% "")
    df <- df[nzchar(trimws(df$smiles)), , drop = FALSE]
    df
  })

  observeEvent(input$run, {
    if (length(problems)) {
      showNotification("Setup is incomplete — run setup.sh first.", type = "error"); return()
    }
    df <- current_input()
    if (!nrow(df)) { showNotification("No SMILES to run.", type = "warning"); return() }
    if (nrow(df) > 500) { showNotification("Please submit at most 500 molecules per batch.", type = "error"); return() }

    rv$results <- NULL; rv$summary <- NULL
    n <- nrow(df)
    prog <- shiny::Progress$new(session, min = 0, max = n)
    prog$set(message = sprintf("Predicting %d molecule(s)…", n), value = 0)
    done_env <- new.env(); done_env$k <- 0
    cb <- function(done, total, name) {
      done_env$k <- done_env$k + 1
      prog$set(value = done_env$k,
               detail = sprintf("%d of %d: %s", done_env$k, total, name))
    }
    on.exit(prog$close(), add = TRUE)

    out <- tryCatch(
      bt_run_batch(df, app_dir = app_dir, btType = input$btType,
                   steps = input$steps, cyp_mode = input$cyp_mode %||% "1",
                   timeout = input$timeout, offline = isTRUE(input$offline),
                   workers = input$workers, progress = cb),
      error = function(e) { showNotification(paste("Run failed:", conditionMessage(e)),
                                             type = "error", duration = NULL); NULL })
    if (is.null(out)) return()
    rv$results <- out$results
    rv$summary <- out$summary
    updateTabsetPanel(session, "tabs", selected = "Metabolites")
  })

  output$has_results <- reactive(!is.null(rv$summary))
  outputOptions(output, "has_results", suspendWhenHidden = FALSE)

  output$status_tiles <- renderUI({
    s <- rv$summary; if (is.null(s)) return(NULL)
    tcount <- function(st) sum(s$Status == st)
    tile <- function(lab, val, cls)
      sprintf('<span class="bt-tile"><b class="%s" style="padding:0 4px;border-radius:6px;">%d</b> %s</span>',
              cls, val, lab)
    div(class = "bt-tiles", style = "margin-bottom:8px;",
        HTML(paste0(
          tile("compounds", nrow(s), "bt-none"),
          tile("predicted", tcount("predicted"), "bt-ok"),
          tile("metabolites", sum(s$Metabolites), "bt-ok"),
          tile("no metabolites", tcount("no metabolites"), "bt-none"),
          tile("parse errors", tcount("parse error"), "bt-parse"),
          tile("errors/timeouts", tcount("error") + tcount("timeout"), "bt-err")
        )))
  })

  output$result_count <- renderText({
    r <- rv$results
    if (is.null(r) || !nrow(r)) "No metabolites." else sprintf("%d metabolite rows", nrow(r))
  })

  # metabolites table (with optional structure thumbnails) -------------------
  output$results_tbl <- DT::renderDataTable({
    r <- rv$results
    if (is.null(r) || !nrow(r)) return(NULL)
    show_cols <- intersect(c("Query", "SMILES", "Molecular formula", "Major Isotope Mass",
                             "ALogP", "Reaction", "Enzyme(s)", "Biosystem",
                             "Metabolite ID", "Precursor SMILES", "InChIKey",
                             "PUBCHEM_CID", "Synonyms"), names(r))
    d <- r[, show_cols, drop = FALSE]
    if ("Enzyme(s)" %in% names(d)) d[["Enzyme(s)"]] <- gsub("\\n", ", ", d[["Enzyme(s)"]])
    if (isTRUE(input$show_struct) && "SMILES" %in% names(d)) {
      struct <- sprintf('<img class="bt-smiles" data-smiles="%s" class="bt-struct"/>',
                        htmltools::htmlEscape(d$SMILES, attribute = TRUE))
      d <- cbind(Structure = struct, d)
    }
    DT::datatable(
      d, escape = FALSE, rownames = FALSE, selection = "none",
      extensions = "Buttons",
      options = list(
        pageLength = 10, scrollX = TRUE, autoWidth = FALSE,
        dom = "lftip",
        drawCallback = DT::JS("function(){ if(window.btDrawStructures) window.btDrawStructures(); }")
      ))
  }, server = TRUE)

  # per-compound summary table ----------------------------------------------
  output$summary_tbl <- DT::renderDataTable({
    s <- rv$summary; if (is.null(s)) return(NULL)
    disp <- s
    disp$Status <- vapply(s$Status, function(st)
      badge(st, if (!is.na(status_class[st])) status_class[st] else "bt-none"),
      character(1))
    DT::datatable(disp, escape = FALSE, rownames = FALSE, selection = "none",
                  options = list(pageLength = 15, scrollX = TRUE, dom = "ftip"))
  }, server = FALSE)

  output$dl_results <- downloadHandler(
    filename = function() sprintf("biotransformer_metabolites_%s.csv", Sys.Date()),
    content = function(file) {
      r <- rv$results
      if (is.null(r) || !nrow(r)) r <- data.frame(message = "no metabolites predicted")
      utils::write.csv(r, file, row.names = FALSE, na = "")
    })

  output$dl_summary <- downloadHandler(
    filename = function() sprintf("biotransformer_summary_%s.csv", Sys.Date()),
    content = function(file) utils::write.csv(rv$summary, file, row.names = FALSE, na = ""))
}

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

shinyApp(ui, server)
