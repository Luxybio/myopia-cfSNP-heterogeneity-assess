# Replace <LOCAL_DATA_ROOT>/<SERVER_DATA_ROOT> with your own paths.
# Compute WT/MT substructure heterogeneity from selected centroids (six classes x count/length/position).

rm(list = ls())

suppressPackageStartupMessages({
  library(RNAsmc)   # getSubStr()
  library(RRNA)     # dependency of RNAsmc
})


SEL_ROOT <- "<LOCAL_DATA_ROOT>/Sfold-seq-build/sfold-selectstructrue-4/selectstructure"
CT_ROOT  <- "<LOCAL_DATA_ROOT>/Sfold-seq-build/sfold-selectstructrue-4"
OUT_DIR  <- "<LOCAL_DATA_ROOT>/substructure-compare-v2-workbuddy-fixed"
LEGACY_DIS <- "<LOCAL_DATA_ROOT>/substructure-compare/each-dis.txt"

## Six substructure classes: output directory name -> key in object returned by RNAsmc::getSubStr()
SUBTYPES  <- c("bulgeLoop", "externalLoop", "hairpinLoop",
               "internalLoop", "multipbranch", "stemLoop")
CLASS_KEY <- c(bulgeLoop        = "bulgeLoop",
               externalLoop     = "externalLoop",
               hairpinLoop      = "hairpinLoop",
               internalLoop     = "internalLoop",
               multipbranch     = "multi-branchLoop",
               stemLoop         = "stem")

##   "membership" (default, follows the reference comment's intent and matches the third/base block)
##        index = union(WT lengths, MT lengths); fill the length value if present, otherwise 0
##   "byorder" (alternative)
##        pair instances one-to-one after sorting by start; pad the tail of the longer side with 0
LEN_PAIR <- "membership"

WRITE_DETAIL <- TRUE   # Whether to write class-contribution / subunit-detail


LOG <- character(0)
log_msg <- function(...) LOG <<- c(LOG, paste0(...))

dir_ok <- function(p) dir.create(p, recursive = TRUE, showWarnings = FALSE)

read_select <- function(path) {
  d <- read.delim(path, sep = "\t", header = TRUE, quote = "",
                  stringsAsFactors = FALSE, check.names = FALSE,
                  colClasses = "character")
  names(d)[1] <- "centroid"
  for (cc in c("matched_pairs", "sum_prob", "pairs", "dG", "selected")) {
    if (cc %in% names(d)) d[[cc]] <- suppressWarnings(as.numeric(d[[cc]]))
  }
  if (!"sum_prob"      %in% names(d)) d$sum_prob      <- NA_real_
  if (!"matched_pairs" %in% names(d)) d$matched_pairs <- NA_real_
  if (!"dG"            %in% names(d)) d$dG            <- NA_real_
  if (!"selected"      %in% names(d)) d$selected      <- 0
  d$centroid <- trimws(d$centroid)
  d
}

## centroid name("c01.ccentroid.ct") -> matching .ct file name fragment("c01-ccentroid.ct")
cen_to_ct <- function(cen) {
  c1 <- sub("\\.ct$", "", cen)                    # c01.ccentroid / ecentroid
  c2 <- sub("\\.ccentroid$", "-ccentroid", c1)    # c01-ccentroid
  paste0(c2, ".ct")
}

ct_read <- function(path) {
  lines <- readLines(path, warn = FALSE)
  seqlen <- NA_real_
  if (length(lines) >= 2) {
    seqlen <- suppressWarnings(
      as.numeric(strsplit(trimws(lines[2]), "\\s+")[[1]][2]))
  }
  ct <- tryCatch(
    read.table(path, skip = 2, header = FALSE, stringsAsFactors = FALSE,
               blank.lines.skip = TRUE),
    error = function(e) NULL)
  if (is.null(ct) || nrow(ct) == 0) return(NULL)
  if (ncol(ct) < 6) {
    ct <- tryCatch(read.table(path, skip = 2, header = FALSE, fill = TRUE,
                              stringsAsFactors = FALSE, blank.lines.skip = TRUE),
                   error = function(e) NULL)
    if (is.null(ct) || ncol(ct) < 6) return(NULL)
  }
  if (!is.finite(seqlen) || seqlen <= 0) seqlen <- nrow(ct)
  list(ct = ct, seqlen = seqlen)
}

## getSubStr() prints "There is no stems" / SFOLD headers to stdout; suppress them
quiet_substr <- function(ct) {
  res <- NULL
  invisible(capture.output(res <- getSubStr(ct)))
  res
}

## Returns a list: each element is the set of base positions (sorted) contained in one substructure unit, sorted by start
inst_of <- function(s, key) {
  v <- s[[key]]
  if (is.null(v) || length(v) == 0) return(list())
  v <- lapply(v, function(x) sort(unique(suppressWarnings(as.integer(x)))))
  v <- v[vapply(v, length, integer(1)) > 0]
  if (length(v) == 0) return(list())
  st <- vapply(v, min, integer(1))
  v[order(st)]
}

## Length block: index = union of length values from both sides (sorted); fill the length value if present, otherwise 0
pad_len_membership <- function(lw, lm) {
  idx <- sort(unique(c(lw, lm)))
  if (length(idx) == 0) return(matrix(0, 0, 2))
  cbind(ifelse(idx %in% lw, idx, 0),
        ifelse(idx %in% lm, idx, 0))
}

## Length block alternative: pair by start order one-to-one; pad the tail of the longer side with 0
pad_len_byorder <- function(lw, lm) {
  n <- max(length(lw), length(lm))
  if (n == 0) return(matrix(0, 0, 2))
  cbind(c(lw, rep(0, n - length(lw))),
        c(lm, rep(0, n - length(lm))))
}

## Position block: index = union of base positions from both sides (sorted); fill the position if present, otherwise 0
## Same as the third step of the reference program (the original used union(), whose order is not sorted,
## but both columns share the same index, so sorting only reorders rows and does not change Euclidean distance)
pad_pos <- function(pw, pm) {
  idx <- sort(union(pw, pm))
  if (length(idx) == 0) return(matrix(0, 0, 2))
  cbind(ifelse(idx %in% pw, idx, 0),
        ifelse(idx %in% pm, idx, 0))
}

stack2 <- function(L) {
  L <- L[!vapply(L, function(x) is.null(x) || nrow(x) == 0, TRUE)]
  if (length(L) == 0) return(matrix(0, 0, 2))
  do.call(rbind, L)
}

eu  <- function(a, b) if (length(a) == 0) 0 else sqrt(sum((a - b)^2))
ruz <- function(a, b) {                       # Ruzicka/Tanimoto, for reference only
  d <- sum(pmax(a, b))
  if (d == 0) return(0)
  1 - sum(pmin(a, b)) / d
}
shannon_H <- function(p) {
  p <- p[is.finite(p) & p > 0]
  if (length(p) == 0) return(0)
  p <- p / sum(p)
  -sum(p * log(p))
}


dir_ok(OUT_DIR)
for (d in c("Num", "Len-subunit-hebing", "Base-subunit-hebing",
            "Num-standard-hebing", "Len-standard-hebing", "Base-standard-hebing")) {
  dir_ok(file.path(OUT_DIR, d))
}
for (st in SUBTYPES) {
  dir_ok(file.path(OUT_DIR, "Len",  st))
  dir_ok(file.path(OUT_DIR, "Base", st))
}


sel_files <- sort(list.files(file.path(SEL_ROOT, "WT"),
                             pattern = "-selectstructure\\.txt$"))
if (length(sel_files) == 0) stop("no selectstructure file under WT dir: ", SEL_ROOT)

pairs <- list()

for (fn in sel_files) {

  key_w <- sub("-selectstructure\\.txt$", "", fn)
  key_m <- sub("^WT_", "MT_", key_w)
  f_m   <- file.path(SEL_ROOT, "MT", paste0(key_m, "-selectstructure.txt"))
  if (!file.exists(f_m)) {
    log_msg("[SKIP] ", key_w, " : MT selectstructure file missing"); next
  }

  dw <- read_select(file.path(SEL_ROOT, "WT", fn))
  dm <- read_select(f_m)

  ## Representative centroid: selected==1; if no marker, take the one with the largest sum_prob
  pick <- function(d, type, key) {
    i <- which(!is.na(d$selected) & d$selected == 1)
    if (length(i) == 0) {
      i <- which.max(ifelse(is.na(d$sum_prob), -Inf, d$sum_prob))
      log_msg("[WARN] ", type, " ", key, " : no selected=1, fallback to max sum_prob: ",
              d$centroid[i])
    }
    if (length(i) > 1) {
      j <- i[which.max(ifelse(is.na(d$sum_prob[i]), -Inf, d$sum_prob[i]))]
      log_msg("[WARN] ", type, " ", key, " : multiple selected=1 (", length(i),
              " rows), keep max sum_prob: ", d$centroid[j])
      i <- j
    }
    i
  }
  iw <- pick(dw, "WT", key_w)
  im <- pick(dm, "MT", key_m)

  ## .ct of that centroid; if missing, fall back in descending sum_prob order
  find_ct <- function(d, i0, type, key) {
    ord <- order(-ifelse(is.na(d$sum_prob), -Inf, d$sum_prob))
    for (j in c(i0, ord[ord != i0])) {
      p <- file.path(CT_ROOT, type, paste0(key, "_", cen_to_ct(d$centroid[j])))
      if (file.exists(p)) {
        if (j != i0) {
          log_msg("[WARN] ", type, " ", key,
                  " : .ct of selected centroid missing, fallback to ", d$centroid[j])
        }
        return(list(path = p, row = j))
      }
    }
    NULL
  }
  cw <- find_ct(dw, iw, "WT", key_w)
  cm <- find_ct(dm, im, "MT", key_m)
  if (is.null(cw) || is.null(cm)) {
    log_msg("[SKIP] ", key_w, " : no .ct found for WT or MT"); next
  }
  iw <- cw$row; im <- cm$row

  pw <- ct_read(cw$path); pm <- ct_read(cm$path)
  if (is.null(pw) || is.null(pm)) {
    log_msg("[SKIP] ", key_w, " : .ct parse failed"); next
  }

  sw <- tryCatch(quiet_substr(pw$ct), error = function(e) NULL)
  sm <- tryCatch(quiet_substr(pm$ct), error = function(e) NULL)
  if (is.null(sw) || is.null(sm)) {
    log_msg("[SKIP] ", key_w, " : getSubStr failed"); next
  }

  numw <- numm <- integer(0)
  lenb <- list(); baseb <- list(); instb <- list()
  for (st in SUBTYPES) {
    k  <- CLASS_KEY[[st]]
    iw_inst <- inst_of(sw, k)
    im_inst <- inst_of(sm, k)
    lw <- vapply(iw_inst, length, integer(1))    # lengths of each instance
    lm <- vapply(im_inst, length, integer(1))
    pwn <- unlist(iw_inst, use.names = FALSE)    # all base positions in this class
    pmn <- unlist(im_inst, use.names = FALSE)

    numw <- c(numw, length(iw_inst))
    numm <- c(numm, length(im_inst))
    lenb[[st]]  <- if (LEN_PAIR == "byorder") pad_len_byorder(lw, lm)
                   else pad_len_membership(lw, lm)
    baseb[[st]] <- pad_pos(pwn, pmn)
    instb[[st]] <- list(w = iw_inst, m = im_inst)
  }

  ## Ensemble-level information (from selectstructure, not used in distance calculation)
  meta <- list(
    cen_w = dw$centroid[iw], cen_m = dm$centroid[im],
    n_w   = nrow(dw), n_m = nrow(dm),
    p_sel_w = dw$sum_prob[iw], p_sel_m = dm$sum_prob[im],
    p_tot_w = sum(dw$sum_prob, na.rm = TRUE),
    p_tot_m = sum(dm$sum_prob, na.rm = TRUE),
    dG_w = dw$dG[iw], dG_m = dm$dG[im],
    H_w = shannon_H(dw$sum_prob), H_m = shannon_H(dm$sum_prob),
    pair_w = dw$matched_pairs[iw], pair_m = dm$matched_pairs[im]
  )

  pairs[[key_w]] <- list(key = key_w, key_m = key_m,
                         seqlen_w = pw$seqlen, seqlen_m = pm$seqlen,
                         numw = numw, numm = numm,
                         lenb = lenb, baseb = baseb, inst = instb,
                         meta = meta)
}

if (length(pairs) == 0) stop("no gene-SNP pair succeeded, check SEL_ROOT / CT_ROOT")
log_msg("[INFO] paired ", length(pairs), " / ", length(sel_files), " gene-SNP items")


dis_tab  <- data.frame()
cls_tab  <- data.frame()
det_tab  <- data.frame()

for (key in names(pairs)) {

  p <- pairs[[key]]

  num_raw <- cbind(p$numw, p$numm)
  storage.mode(num_raw) <- "double"
  ns <- max(num_raw)
  if (!is.finite(ns) || ns <= 0) ns <- 1
  num_std <- num_raw / ns

  ## ---- 4.2 Length: after merging six classes, divide by the per-sample merged matrix maximum ----
  len_raw <- stack2(p$lenb)
  ls <- max(len_raw)
  if (!is.finite(ls) || ls <= 0) ls <- 1
  len_std <- len_raw / ls

  base_raw <- stack2(p$baseb)
  base_std <- cbind(base_raw[, 1] / p$seqlen_w,
                    base_raw[, 2] / p$seqlen_m)

  write.table(num_raw,  file.path(OUT_DIR, "Num", paste0(key, ".txt")),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  for (st in SUBTYPES) {
    write.table(p$lenb[[st]],  file.path(OUT_DIR, "Len",  st, paste0(key, ".txt")),
                sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    write.table(p$baseb[[st]], file.path(OUT_DIR, "Base", st, paste0(key, ".txt")),
                sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  write.table(num_std,  file.path(OUT_DIR, "Num-standard-hebing",  paste0(key, ".txt")),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(len_raw,  file.path(OUT_DIR, "Len-subunit-hebing",   paste0(key, ".txt")),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(base_raw, file.path(OUT_DIR, "Base-subunit-hebing",  paste0(key, ".txt")),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(len_std,  file.path(OUT_DIR, "Len-standard-hebing",  paste0(key, ".txt")),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(base_std, file.path(OUT_DIR, "Base-standard-hebing", paste0(key, ".txt")),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

  ## ---- 4.5 Distance: rbind the three blocks, then WT column vs MT column, ordinary Euclidean distance ----
  allcon <- rbind(num_std, len_std, base_std)
  d_eu   <- eu(allcon[, 1], allcon[, 2])
  d_num  <- eu(num_std[, 1],  num_std[, 2])
  d_len  <- eu(len_std[, 1],  len_std[, 2])
  d_base <- eu(base_std[, 1], base_std[, 2])
  d_ruz  <- ruz(allcon[, 1], allcon[, 2])

  m <- p$meta
  dis_tab <- rbind(dis_tab, data.frame(
    id = key, id_m = p$key_m,
    euclidean = d_eu,
    d_num = d_num, d_len = d_len, d_base = d_base,
    ruzicka = d_ruz,
    dim_num = nrow(num_std), dim_len = nrow(len_std), dim_base = nrow(base_std),
    seqlen_w = p$seqlen_w, seqlen_m = p$seqlen_m,
    cen_w = m$cen_w, cen_m = m$cen_m, n_cen_w = m$n_w, n_cen_m = m$n_m,
    prob_sel_w = m$p_sel_w, prob_sel_m = m$p_sel_m,
    prob_share_w = if (is.finite(m$p_tot_w) && m$p_tot_w > 0) m$p_sel_w / m$p_tot_w else NA,
    prob_share_m = if (is.finite(m$p_tot_m) && m$p_tot_m > 0) m$p_sel_m / m$p_tot_m else NA,
    H_w = m$H_w, H_m = m$H_m,
    expH_w = exp(m$H_w), expH_m = exp(m$H_m),
    dG_w = m$dG_w, dG_m = m$dG_m,
    dG_delta = m$dG_m - m$dG_w,
    paired_n_w = m$pair_w, paired_n_m = m$pair_m,
    stringsAsFactors = FALSE))

  if (WRITE_DETAIL) {
    for (st in SUBTYPES) {
      lb <- p$lenb[[st]]; bb <- p$baseb[[st]]
      ls_st <- lb / ls
      bs_st <- cbind(bb[, 1] / p$seqlen_w, bb[, 2] / p$seqlen_m)
      cls_tab <- rbind(cls_tab, data.frame(
        id = key, class = st,
        n_w = p$numw[match(st, SUBTYPES)], n_m = p$numm[match(st, SUBTYPES)],
        len_sum_w = if (nrow(lb) > 0) sum(lb[, 1]) else 0,
        len_sum_m = if (nrow(lb) > 0) sum(lb[, 2]) else 0,
        d_len  = eu(ls_st[, 1], ls_st[, 2]),
        d_base = eu(bs_st[, 1], bs_st[, 2]),
        stringsAsFactors = FALSE))
      for (side in c("w", "m")) {
        ins <- p$inst[[st]][[side]]
        if (length(ins) == 0) next
        det_tab <- rbind(det_tab, data.frame(
          id = key, type = if (side == "w") "WT" else "MT", class = st,
          index = seq_along(ins),
          start = vapply(ins, min, integer(1)),
          end   = vapply(ins, max, integer(1)),
          length = vapply(ins, length, integer(1)),
          positions = vapply(ins, function(x) paste(x, collapse = ","), ""),
          stringsAsFactors = FALSE))
      }
    }
  }

  log_msg(sprintf("[OK] %-62s eu=%.4f (num %.4f / len %.4f / base %.4f) dims %d/%d/%d",
                  key, d_eu, d_num, d_len, d_base,
                  nrow(num_std), nrow(len_std), nrow(base_std)))
}


dis_tab <- dis_tab[order(-dis_tab$euclidean), ]
dis_tab <- cbind(rank = seq_len(nrow(dis_tab)), dis_tab)

write.table(dis_tab, file.path(OUT_DIR, "each-dis.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

if (WRITE_DETAIL) {
  write.table(cls_tab, file.path(OUT_DIR, "class-contribution.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(det_tab, file.path(OUT_DIR, "subunit-detail.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

write.table(dis_tab[, c("id", "id_m", "n_cen_w", "n_cen_m", "H_w", "H_m",
                        "expH_w", "expH_m", "prob_share_w", "prob_share_m",
                        "dG_w", "dG_m", "dG_delta", "euclidean")],
            file.path(OUT_DIR, "ensemble-heterogeneity.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

if (file.exists(LEGACY_DIS)) {
  lg <- tryCatch(read.delim(LEGACY_DIS, sep = "\t", stringsAsFactors = FALSE),
                 error = function(e) NULL)
  if (!is.null(lg) && "id" %in% names(lg)) {
    lg$id <- sub("_c[0-9]+-ccentroid$", "", lg$id)
    lg <- lg[!duplicated(lg$id), intersect(c("id", "euclidean", "d_num", "d_len",
                                             "d_base"), names(lg))]
    names(lg)[-1] <- paste0(names(lg)[-1], "_legacy")
    cur <- dis_tab[, c("id", "euclidean", "d_num", "d_len", "d_base")]
    names(cur)[-1] <- paste0(names(cur)[-1], "_V3")
    cmpdf <- merge(lg, cur, by = "id", all = FALSE)
    if (nrow(cmpdf) > 0) {
      cmpdf$rank_legacy <- rank(-cmpdf$euclidean_legacy, ties.method = "min")
      cmpdf$rank_V3     <- rank(-cmpdf$euclidean_V3,     ties.method = "min")
      cmpdf$rank_shift  <- cmpdf$rank_legacy - cmpdf$rank_V3
      cmpdf <- cmpdf[order(-cmpdf$euclidean_V3), ]
      write.table(cmpdf, file.path(OUT_DIR, "each-dis-vs-legacy.txt"),
                  sep = "\t", quote = FALSE, row.names = FALSE)
      r <- suppressWarnings(cor(cmpdf$euclidean_legacy, cmpdf$euclidean_V3,
                                method = "spearman"))
      log_msg(sprintf("[INFO] legacy comparison on %d items, Spearman rho = %.3f",
                      nrow(cmpdf), r))
    }
  }
}

con <- file(file.path(OUT_DIR, "substructure-compare-V3.log.txt"), "w")
writeLines(LOG, con); close(con)


cat(sprintf("\nDONE: %d gene-SNP items -> %s\n", nrow(dis_tab), OUT_DIR))
cat(sprintf("params: LEN_PAIR=%s | normalisation per-sample (num by max count, len by max length, base by seqlen)\n",
            LEN_PAIR))
cat("\nSummary of euclidean:\n"); print(summary(dis_tab$euclidean))
cat("\nTop 10 by euclidean distance:\n")
print(head(dis_tab[, c("id", "euclidean", "d_num", "d_len", "d_base",
                       "dim_num", "dim_len", "dim_base")], 10), row.names = FALSE)
cat(sprintf("\nlog: %s\n", file.path(OUT_DIR, "substructure-compare-V3.log.txt")))


