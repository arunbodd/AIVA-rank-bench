#!/usr/bin/env Rscript
# _clinprior_worker.R  –  DO NOT RUN DIRECTLY; called by 08_run_clinprior.py
#
# Processes a batch of cases with ClinPrior in phenotype-only mode:
#   proteinScore()  →  MatrixPropagation()  →  combined rank
#
# Args:  <input_jsonl>  <output_jsonl>  <top_n>

suppressPackageStartupMessages({
  library(ClinPrior)
  library(jsonlite)
})

args       <- commandArgs(trailingOnly = TRUE)
input_path  <- args[1]
output_path <- args[2]
top_n       <- as.integer(args[3])

# Mitochondrial gene symbol mapping: HGNC MT- prefix → ClinPrior symbol
MT_GENE_MAP <- c(
  "MT-ATP6" = "ATP6",  "MT-ATP8" = "ATP8",
  "MT-CO1"  = "COX1",  "MT-CO2"  = "COX2",  "MT-CO3"  = "COX3",
  "MT-CYB"  = "CYTB",
  "MT-ND1"  = "ND1",   "MT-ND2"  = "ND2",   "MT-ND3"  = "ND3",
  "MT-ND4"  = "ND4",   "MT-ND4L" = "ND4L",  "MT-ND5"  = "ND5",  "MT-ND6" = "ND6",
  "MT-TA"   = "TRNA",  "MT-TC"   = "TRNC",  "MT-TD"   = "TRND",
  "MT-TE"   = "TRNE",  "MT-TF"   = "TRNF",  "MT-TG"   = "TRNG",
  "MT-TH"   = "TRNH",  "MT-TI"   = "TRNI",  "MT-TK"   = "TRNK",
  "MT-TL1"  = "TRNL1", "MT-TL2"  = "TRNL2", "MT-TM"   = "TRNM",
  "MT-TN"   = "TRNN",  "MT-TP"   = "TRNP",  "MT-TQ"   = "TRNQ",
  "MT-TR"   = "TRNR",  "MT-TS1"  = "TRNS1", "MT-TS2"  = "TRNS2",
  "MT-TT"   = "TRNT",  "MT-TV"   = "TRNV",  "MT-TW"   = "TRNW",
  "MT-TY"   = "TRNY"
)

resolve_gene <- function(symbol) {
  mapped <- MT_GENE_MAP[symbol]
  if (!is.na(mapped)) mapped else symbol
}

cases <- lapply(readLines(input_path), fromJSON)
total <- length(cases)

con_out <- file(output_path, open = "wt")

for (i in seq_along(cases)) {
  case      <- cases[[i]]
  pid       <- case$patient_id
  hpo       <- unlist(case$observed_hpo)
  true_gene       <- as.character(case$true_gene_symbol)
  true_gene_query <- resolve_gene(true_gene)   # ClinPrior network symbol

  cat(sprintf("[%d/%d] pid=%s  hpo_n=%d\n", i, total, pid, length(hpo)),
      file = stderr())

  result_base <- list(
    patient_id           = pid,
    disease_id           = case$disease_id,
    broad_category       = case$broad_category,
    true_gene_symbol     = true_gene,
    clinprior_rank       = NULL,
    clinprior_rank_funct = NULL,
    clinprior_rank_phys  = NULL,
    clinprior_n_genes    = NA_integer_,
    clinprior_top_genes  = list()
  )

  # ── Step 1: HPO → gene scores ───────────────────────────────────────────────
  score_result <- tryCatch(
    proteinScore(hpo, verbose = FALSE),
    error = function(e) {
      cat(sprintf("  [ERROR] proteinScore: %s\n", conditionMessage(e)),
          file = stderr())
      NULL
    }
  )
  if (is.null(score_result)) {
    result_base$error <- "proteinScore failed"
    writeLines(toJSON(result_base, auto_unbox = TRUE, null = "null"), con_out)
    flush(con_out)
    next
  }

  # ── Step 2: Propagate through gene networks ─────────────────────────────────
  gps <- tryCatch(
    MatrixPropagation(score_result$scores, alpha = 0.2, verbose = FALSE),
    error = function(e) {
      cat(sprintf("  [ERROR] MatrixPropagation: %s\n", conditionMessage(e)),
          file = stderr())
      NULL
    }
  )
  if (is.null(gps)) {
    result_base$error <- "MatrixPropagation failed"
    writeLines(toJSON(result_base, auto_unbox = TRUE, null = "null"), con_out)
    flush(con_out)
    next
  }

  # ── Step 3: Merge functional + physical scores, rank genes ──────────────────
  df_funct <- data.frame(
    gene        = as.character(gps$Symbol_Funct),
    score_funct = as.numeric(gps$PriorFunct),
    stringsAsFactors = FALSE
  )
  df_phys <- data.frame(
    gene       = as.character(gps$Symbol_Phys),
    score_phys = as.numeric(gps$PriorPhys),
    stringsAsFactors = FALSE
  )

  merged <- merge(df_funct, df_phys, by = "gene", all = TRUE)
  merged$score_funct[is.na(merged$score_funct)] <- 0
  merged$score_phys [is.na(merged$score_phys)]  <- 0
  merged$combined   <- (merged$score_funct + merged$score_phys) / 2
  merged            <- merged[order(-merged$combined), ]

  n_genes       <- nrow(merged)
  rank_combined <- match(true_gene_query, merged$gene)   # NA if absent
  rank_funct    <- match(true_gene_query, df_funct$gene)
  rank_phys     <- match(true_gene_query, df_phys$gene)

  # Top-N gene entries for inspection
  top_k     <- min(top_n, nrow(merged))
  top_genes <- lapply(seq_len(top_k), function(j) list(
    gene           = merged$gene[j],
    score_combined = round(merged$combined[j],    6),
    score_funct    = round(merged$score_funct[j], 6),
    score_phys     = round(merged$score_phys[j],  6)
  ))

  record <- list(
    patient_id           = pid,
    disease_id           = case$disease_id,
    broad_category       = case$broad_category,
    true_gene_symbol     = true_gene,
    clinprior_rank       = if (is.na(rank_combined)) NULL else as.integer(rank_combined),
    clinprior_rank_funct = if (is.na(rank_funct))    NULL else as.integer(rank_funct),
    clinprior_rank_phys  = if (is.na(rank_phys))     NULL else as.integer(rank_phys),
    clinprior_n_genes    = n_genes,
    clinprior_top_genes  = top_genes
  )

  writeLines(toJSON(record, auto_unbox = TRUE, null = "null"), con_out)
  flush(con_out)
}

close(con_out)
cat(sprintf("[ClinPrior worker] done: %d cases written to %s\n", total, output_path),
    file = stderr())
