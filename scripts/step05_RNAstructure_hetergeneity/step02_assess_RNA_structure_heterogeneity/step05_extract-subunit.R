# Replace <LOCAL_DATA_ROOT>/<SERVER_DATA_ROOT> with your own paths.
# Decompose each WT/MT centroid .ct into six substructure classes and write per-class/per-sample files plus merged summaries.

rm(list = ls())

suppressPackageStartupMessages({
  library(RNAsmc)   # getSubStr()
  library(RRNA)     # required by RNAsmc
})

IN_ROOT  <- "<LOCAL_DATA_ROOT>/Sfold-seq-build/sfold-selectstructrue-4"
OUT_ROOT <- "<LOCAL_DATA_ROOT>/substructure-compute"
TYPES    <- c("WT", "MT")

## class keys in RNAsmc's return object -> output sub-directory name
CLASS_MAP <- c(
  bulgeLoop          = "bulgeLoop",
  externalLoop       = "externalLoop",
  hairpinLoop        = "hairpinLoop",
  internalLoop       = "internalLoop",
  `multip-branchLoop`= "multipbranch",
  stem               = "stemLoop"
)
## label written in column 1 of each class file
CLASS_LABEL <- c(
  bulgeLoop          = "bulgeLoop",
  externalLoop       = "externalLoop",
  hairpinLoop        = "hairpinLoop",
  internalLoop       = "internalLoop",
  `multip-branchLoop`= "multipbranchLoop",
  stem               = "stem"
)
## order of numArr returned by getSubStr()
NUM_ORDER <- c("bulgeLoop", "externalLoop", "hairpinLoop",
               "internalLoop", "multip-branchLoop", "stem")

LOG <- character(0)
log_msg <- function(...) LOG <<- c(LOG, paste0(...))

read_ct <- function(f) {
  lines <- readLines(f, warn = FALSE)
  hdr_len <- NA_integer_
  if (length(lines) >= 2 && grepl("^Length", lines[2])) {
    hdr_len <- suppressWarnings(
      as.integer(sub("^Length:\\s*([0-9]+).*$", "\\1", lines[2]))
    )
  }
  ct <- read.table(f, skip = 2, header = FALSE, stringsAsFactors = FALSE,
                   blank.lines.skip = TRUE)
  list(ct = ct, hdr_len = hdr_len)
}

## getSubStr prints "There is no stems" for stem-free structures: swallow it
quiet_substr <- function(ct) {
  invisible(capture.output(res <- getSubStr(ct)))
  res
}

## write one class file (same 5-field layout as the original script)
write_class_file <- function(vals, label, id, path) {
  vals <- vals[vapply(vals, length, integer(1)) > 0]
  if (length(vals) == 0) {
    line <- paste(label, "0", "0", "0", id, sep = "\t")
  } else {
    pos <- paste(unlist(vals), collapse = ",")
    len <- paste(vapply(vals, length, integer(1)), collapse = ",")
    tot <- length(unlist(vals))
    line <- paste(label, pos, len, tot, id, sep = "\t")
  }
  writeLines(line, path, useBytes = TRUE)
}

dir.create(OUT_ROOT, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_ROOT, "summary"), recursive = TRUE, showWarnings = FALSE)

all_num  <- list()   # per type: data.frame of class counts
all_long <- list()   # per type: data.frame of one row per substructure

for (type in TYPES) {

  in_dir   <- file.path(IN_ROOT, type)
  feat_dir <- file.path(OUT_ROOT, type, "sig-substructure-feacture")
  num_dir  <- file.path(OUT_ROOT, type, "sig-substructure-num")

  if (!dir.exists(in_dir)) {
    log_msg("[FATAL] input dir not found: ", in_dir)
    next
  }
  for (d in CLASS_MAP) dir.create(file.path(feat_dir, d), recursive = TRUE,
                                  showWarnings = FALSE)
  dir.create(num_dir, recursive = TRUE, showWarnings = FALSE)

  files <- list.files(in_dir, pattern = "\\.ct$", full.names = FALSE)
  files <- sort(files)
  cat(sprintf("[%s] %d .ct files in %s\n", type, length(files), in_dir))

  num_rows  <- list()
  long_rows <- list()

  for (fname in files) {
    id <- tools::file_path_sans_ext(fname)          # keeps dots inside the ID
    fpath <- file.path(in_dir, fname)

    parsed <- tryCatch(read_ct(fpath), error = function(e) NULL)
    if (is.null(parsed)) {
      log_msg("[ERROR] ", type, "\t", id, "\tread failed")
      next
    }
    ct <- parsed$ct
    if (ncol(ct) != 6) {
      log_msg("[WARN] ", type, "\t", id, "\texpected 6 columns, got ", ncol(ct))
    }
    if (!is.na(parsed$hdr_len) && nrow(ct) != parsed$hdr_len) {
      log_msg("[WARN] ", type, "\t", id, "\theader length ", parsed$hdr_len,
              " != ", nrow(ct), " data rows")
    }

    s <- tryCatch(quiet_substr(ct), error = function(e) {
      log_msg("[ERROR] ", type, "\t", id, "\tgetSubStr failed: ",
              conditionMessage(e))
      NULL
    })
    if (is.null(s)) next

    for (k in names(CLASS_MAP)) {
      write_class_file(s[[k]], CLASS_LABEL[[k]], id,
                       file.path(feat_dir, CLASS_MAP[[k]], paste0(id, ".txt")))
    }

    counts <- vapply(NUM_ORDER, function(k) length(s[[k]]), integer(1))
    if (!is.null(s$numArr) && !identical(as.integer(counts), as.integer(s$numArr))) {
      log_msg("[WARN] ", type, "\t", id, "\tnumArr mismatch: recomputed ",
              paste(counts, collapse = ","), " vs numArr ",
              paste(s$numArr, collapse = ","))
    }
    line <- paste(c(id, counts), collapse = "\t")
    writeLines(line, file.path(num_dir, paste0(id, ".txt")), useBytes = TRUE)

    num_rows[[id]] <- data.frame(
      ID = id, as.list(counts), check.names = FALSE, stringsAsFactors = FALSE
    )

    for (k in names(CLASS_MAP)) {
      v <- s[[k]]
      if (length(v) == 0) next
      for (j in seq_along(v)) {
        p <- sort(as.integer(v[[j]]))
        p <- p[!is.na(p)]
        if (length(p) == 0) next
        long_rows[[length(long_rows) + 1]] <- data.frame(
          ID         = id,
          class      = CLASS_LABEL[[k]],
          index      = j,
          start      = min(p),
          end        = max(p),
          length     = length(p),
          positions  = paste(p, collapse = ","),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (length(num_rows) > 0) {
    num_df <- do.call(rbind, num_rows)
    rownames(num_df) <- NULL
    names(num_df) <- c("ID", NUM_ORDER)
    write.table(num_df,
                file.path(OUT_ROOT, "summary",
                          paste0(type, "-substructure-num.txt")),
                sep = "\t", quote = FALSE, row.names = FALSE)
    all_num[[type]] <- num_df
  }
  if (length(long_rows) > 0) {
    long_df <- do.call(rbind, long_rows)
    rownames(long_df) <- NULL
    write.table(long_df,
                file.path(OUT_ROOT, "summary",
                          paste0(type, "-substructure-long.txt")),
                sep = "\t", quote = FALSE, row.names = FALSE)
    all_long[[type]] <- long_df
  }

  cat(sprintf("[%s] done: %d samples, %d substructure rows\n",
              type, length(num_rows), length(long_rows)))
}

log_path <- file.path(OUT_ROOT, "substructure-extract.log.txt")
if (length(LOG) == 0) {
  writeLines("[OK] no warnings or errors", log_path, useBytes = TRUE)
} else {
  writeLines(LOG, log_path, useBytes = TRUE)
}
cat("Log:", log_path, "\n")
cat("Output root:", OUT_ROOT, "\n")
