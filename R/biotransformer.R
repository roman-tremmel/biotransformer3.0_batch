# -----------------------------------------------------------------------------
# BioTransformer 3.0 batch engine wrapper
#
# Runs the official BioTransformer 3.0 runnable JAR once per input molecule so
# that:
#   * every input keeps its own row in the results, even when it cannot be
#     parsed or predicted (the batch never aborts on a single bad molecule);
#   * every predicted metabolite is attributed to the exact query it came from.
#
# Only base R + processx are required (no RDKit / no extra chemistry stack).
# -----------------------------------------------------------------------------

BT_OUTPUT_COLS <- c(
  "InChI", "InChIKey", "SMILES", "Synonyms", "PUBCHEM_CID", "Molecular formula",
  "Major Isotope Mass", "ALogP", "Lipinski_Violations",
  "Insecticide_Likeness_Violations", "Post_Em_Herbicide_Likeness_Violations",
  "Metabolite ID", "cdk:Title", "Reaction", "Reaction ID", "Enzyme(s)",
  "Biosystem", "Precursor ID", "Precursor SMILES", "Precursor InChI",
  "Precursor InChIKey", "Precursor ALogP", "Precursor Major Isotope Mass"
)

# Biotransformer types accepted by the JAR (-b). Labels shown in the UI.
BT_TYPES <- c(
  "All human (recommended)"       = "allHuman",
  "CYP450"                        = "cyp450",
  "Phase II"                      = "phaseII",
  "Human gut microbial"           = "hgut",
  "EC-based"                      = "ecbased",
  "Human super-transformer"       = "superbio",
  "Environmental microbial"       = "env"
)

BT_CYP_MODES <- c(
  "1 - CypReact + BioTransformer rules"            = "1",
  "2 - CyProducts only"                            = "2",
  "3 - Combined (CypReact + rules + CyProducts)"   = "3"
)

# ---- locating the runnable JAR --------------------------------------------

#' Find the BioTransformer JAR to run.
#' @param app_dir install directory (holds config.json, btkb/, supportfiles/)
#' @param offline prefer the offline (no online annotation) variant if present
bt_find_jar <- function(app_dir = ".", offline = FALSE) {
  cand <- character(0)
  if (offline) cand <- c(cand, file.path(app_dir, "BioTransformer3.0.offline.jar"))
  cand <- c(cand, file.path(app_dir, "BioTransformer3.0.jar"))
  # fall back to any BioTransformer*.jar (e.g. the raw vendor download)
  cand <- c(cand, sort(Sys.glob(file.path(app_dir, "BioTransformer3.0*.jar")),
                        decreasing = TRUE))
  cand <- cand[file.exists(cand)]
  if (!length(cand)) return(NA_character_)
  cand[[1]]
}

#' Quick sanity check that the install directory is complete.
bt_check_install <- function(app_dir = ".") {
  problems <- character(0)
  if (is.na(bt_find_jar(app_dir)))
    problems <- c(problems, "BioTransformer JAR not found (run setup.sh).")
  for (f in c("config.json", "btkb", "supportfiles"))
    if (!file.exists(file.path(app_dir, f)))
      problems <- c(problems, sprintf("missing '%s' next to the JAR.", f))
  if (Sys.which("java") == "")
    problems <- c(problems, "java not found on PATH.")
  problems
}

# ---- running one molecule --------------------------------------------------

.bt_clean_msg <- function(x) {
  if (!length(x)) return("")
  lines <- unlist(strsplit(paste(x, collapse = "\n"), "\n", fixed = TRUE))
  drop <- grepl(paste0("packet sent successfully|Communications link failure|",
                       "CJCommunicationsException|CommunicationsException|",
                       "^\\s*at (com\\.mysql|java\\.|sun\\.)|^\\s*\\.\\.\\.|",
                       "Caused by:|Picked up JAVA_TOOL_OPTIONS"),
                lines)
  paste(trimws(lines[!drop & nzchar(trimws(lines))]), collapse = " | ")
}

#' Run BioTransformer for a single molecule.
#'
#' @return list(status, n_metabolites, results (data.frame or NULL), message,
#'              seconds)
bt_run_one <- function(smiles, name = NA, id = NA,
                       app_dir = ".",
                       jar = NULL,
                       btType = "allHuman",
                       steps = 1L,
                       cyp_mode = "1",
                       timeout = 300,
                       xmx = NULL) {
  smiles <- trimws(as.character(smiles))
  if (is.null(jar)) jar <- bt_find_jar(app_dir)
  res <- list(status = "error", n_metabolites = 0L, results = NULL,
              message = "", seconds = NA_real_)
  if (is.na(jar) || !nzchar(jar)) { res$message <- "JAR not found"; return(res) }
  if (!nzchar(smiles)) { res$status <- "parse error"; res$message <- "empty SMILES"; return(res) }

  work <- tempfile("bt_")
  dir.create(work, recursive = TRUE, showWarnings = FALSE)
  jhome <- file.path(work, "home")
  dir.create(jhome, showWarnings = FALSE)
  out_csv <- file.path(work, "out.csv")
  on.exit(unlink(work, recursive = TRUE, force = TRUE), add = TRUE)

  args <- c(
    if (!is.null(xmx)) sprintf("-Xmx%s", xmx),
    # JNI-InChI unpacks a native lib to $HOME/.jnati; give it a writable home
    # (fixes "Error creating directory" on Posit Connect / restricted servers)
    sprintf("-Duser.home=%s", jhome),
    sprintf("-Djava.io.tmpdir=%s", work),
    "-jar", jar,
    "-k", "pred",
    "-b", btType,
    "-ismi", smiles,
    "-ocsv", out_csv,
    "-s", as.character(as.integer(steps))
  )
  if (identical(btType, "allHuman") || identical(btType, "cyp450"))
    args <- c(args, "-cm", as.character(cyp_mode))

  t0 <- Sys.time()
  pr <- tryCatch(
    processx::run("java", args, wd = normalizePath(app_dir),
                  error_on_status = FALSE, timeout = timeout,
                  stderr_to_stdout = FALSE),
    error = function(e) structure(list(timeout = inherits(e, "system_command_timeout_error"),
                                       status = NA_integer_,
                                       stdout = "", stderr = conditionMessage(e)),
                                  class = "bt_err"))
  res$seconds <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)

  combined <- paste(pr$stdout, pr$stderr, sep = "\n")
  if (isTRUE(pr$timeout)) {
    res$status <- "timeout"
    res$message <- sprintf("exceeded %ss time limit", timeout)
    return(res)
  }
  if (grepl("InvalidSmilesException|could not parse|unexpected character", combined)) {
    res$status <- "parse error"
    res$message <- "SMILES could not be parsed"
    return(res)
  }

  if (file.exists(out_csv)) {
    df <- tryCatch(utils::read.csv(out_csv, check.names = FALSE,
                                   stringsAsFactors = FALSE, colClasses = "character"),
                   error = function(e) NULL)
    if (!is.null(df)) {
      if (nrow(df) > 0) {
        res$status <- "predicted"
        res$n_metabolites <- nrow(df)
        res$results <- df
      } else {
        res$status <- "no metabolites"
      }
      return(res)
    }
  }

  # No usable output and not a parse error -> genuine failure
  if (!is.na(pr$status) && pr$status == 0L) {
    res$status <- "no metabolites"
  } else {
    res$status <- "error"
    res$message <- .bt_clean_msg(combined)
    if (!nzchar(res$message)) res$message <- sprintf("java exit status %s", pr$status)
  }
  res
}

# ---- running a batch -------------------------------------------------------

#' Run a batch of molecules.
#'
#' @param df data.frame with a `smiles` column and optional `name` column.
#' @param progress function(done, total, name) or NULL.
#' @param workers number of parallel JVMs (fork-based; 1 = sequential).
#' @return list(results = combined metabolite data.frame with Query_* columns,
#'              summary = per-input status data.frame)
bt_run_batch <- function(df, app_dir = ".", btType = "allHuman", steps = 1L,
                         cyp_mode = "1", timeout = 300, offline = FALSE,
                         workers = 1L, xmx = NULL, progress = NULL) {
  stopifnot("smiles" %in% names(df))
  n <- nrow(df)
  if (!"name" %in% names(df)) df$name <- ""
  df$name <- ifelse(is.na(df$name) | !nzchar(trimws(df$name)),
                    sprintf("mol_%d", seq_len(n)), trimws(df$name))
  jar <- bt_find_jar(app_dir, offline = offline)

  one <- function(i) {
    bt_run_one(df$smiles[i], name = df$name[i], id = i, app_dir = app_dir,
               jar = jar, btType = btType, steps = steps,
               cyp_mode = cyp_mode, timeout = timeout, xmx = xmx)
  }
  # Progress is reported from THIS (main) process only -- a Shiny Progress
  # handle cannot be used inside a forked mclapply worker.
  tick <- function(done, name) if (!is.null(progress)) try(progress(done, n, name), silent = TRUE)

  workers <- max(1L, as.integer(workers))
  results <- vector("list", n)
  if (workers > 1L && .Platform$OS.type == "unix" && n > 1L) {
    # run in chunks of `workers`; update progress in the parent after each chunk
    idx <- split(seq_len(n), ceiling(seq_len(n) / workers))
    done <- 0L
    for (chunk in idx) {
      part <- parallel::mclapply(chunk, one, mc.cores = workers, mc.preschedule = FALSE)
      for (p in seq_along(chunk)) {
        r <- part[[p]]
        if (inherits(r, "try-error") || !is.list(r))
          r <- list(status = "error", n_metabolites = 0L, results = NULL,
                    message = paste("worker error:", trimws(as.character(r))), seconds = NA_real_)
        results[[chunk[p]]] <- r
        done <- done + 1L
        tick(done, df$name[chunk[p]])
      }
    }
  } else {
    for (i in seq_len(n)) {
      r <- tryCatch(one(i), error = function(e)
        list(status = "error", n_metabolites = 0L, results = NULL,
             message = conditionMessage(e), seconds = NA_real_))
      results[[i]] <- r
      tick(i, df$name[i])
    }
  }

  # Assemble combined metabolite table + per-input summary
  met_list <- list()
  summ <- data.frame(Query = df$name, SMILES = df$smiles,
                     Status = NA_character_, Metabolites = 0L,
                     Seconds = NA_real_, Note = "", stringsAsFactors = FALSE)
  for (i in seq_len(n)) {
    r <- results[[i]]
    summ$Status[i]      <- r$status
    summ$Metabolites[i] <- r$n_metabolites
    summ$Seconds[i]     <- r$seconds
    summ$Note[i]        <- r$message
    if (!is.null(r$results) && nrow(r$results) > 0) {
      pre <- data.frame(Query = rep(df$name[i], nrow(r$results)),
                        `Query SMILES` = rep(df$smiles[i], nrow(r$results)),
                        check.names = FALSE, stringsAsFactors = FALSE)
      met_list[[length(met_list) + 1L]] <- cbind(pre, r$results)
    }
  }
  results_df <- if (length(met_list))
    do.call(function(...) rbind(..., make.row.names = FALSE), met_list)
  else data.frame()
  list(results = results_df, summary = summ)
}
