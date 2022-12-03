#' Select hits by RSA
#'
#' Selected hits by redundant siRNA activity method. Here is a wrapper function of RSA 1.8 by Yingyao Zhou.
#'
#' @param dta synthetic lethal RNAi screen data
#' @param treatment the treatment condition in EXPERIMENT_MODIFICATION
#' @param control the control condition in EXPERIMENT_MODIFICATION
#' @param normMethod normalization methods. If "PLATE", then values are normalized by plate median, otherwise use the
#' provided control siRNA
#' @param LB Low bound
#' @param UB up bound
#' @param revHits reverse hit picking, default the lower the score the better
#' @param Bonferroni conceptually useful when there are different number of siRNAs per gene, default FALSE
#' @param outputFile output file name
#' @param scoreFile name of the score file to be written under the current folder
#' @return A result file written to the current folder.
#' \itemize{
#'   \item Gene_ID,Well_ID,Score: columns from input spreadsheet
#'   \item LogP: OPI p-value in log10, i.e., -2 means 0.01
#'   \item OPI_Hit: whether the well is a hit, 1 means yes, 0 means no
#'   \item #hitWell: number of hit wells for the gene
#'   \item #totalWell: total number of wells for the gene. If gene A has three wells w1, w2 and w3, and w1 and w2 are hits, #totalWell should be 3, #hitWell should be 2, w1 and w2 should have OPI_Hit set as 1 and w3 should have OPI_Hit set as 0.
#'   \item OPI_Rank: ranking column to sort all wells for hit picking
#'   \item Cutoff_Rank: ranking column to sort all wells based on Score in the simple activity-based method
#' }
#' Note: a rank value of 999999 means the well is not a hit
#' @references
#' Koenig, R. et al. A probability-based approach for the analysis of large-scale RNAi screens. Nat Methods 4, 847-849 (2007).
#' @examples
#' res <- rsaHits(example_dt, treatment = "treatment", control = "control",
#'                normMethod = "PLATE", LB = 0.2, UB = 0.8, revHits = FALSE,
#'                Bonferroni = FALSE, outputFile = "RSAhits.csv")
#' @export
rsaHits <- function(dta,
    treatment,
    control,
    normMethod = "PLATE",
    LB,
    UB,
    revHits    = FALSE,
    Bonferroni = FALSE,
    outputFile = "RSAhits.csv",
    scoreFile  = "RSA_score.csv") {
  tem.1 <- sapply(unique(as.character(dta$MASTER_PLATE)), .ff_rsaRatio, dta, treatment, control, normMethod = normMethod, simplify = FALSE)
  rsaScore <- data.frame(do.call(rbind, lapply(names(tem.1), function(x) tem.1[[x]])))
  write.table(rsaScore, scoreFile, sep = "\t", quote = FALSE, row.names = FALSE)

  opts <- list(LB = LB, UB = UB, outputFile = outputFile, reverse = revHits, bonferroni = Bonferroni)
  rsaRes <- OPI(rsaScore$Gene_ID, rsaScore$Score, opts, rsaScore)
  message("File written to:", outputFile)
  write.table(rsaRes, outputFile, sep = "\t", quote = FALSE, row.names = FALSE)

  message("Summary of RSA:")
  message(paste("Total #Genes = ", length(unique(rsaRes$Gene_ID)), collapse = ""));
  message(paste("Total #Wells = ", nrow(rsaRes), collapse = ""));
  message(paste("Total #Hit Genes = ", length(unique(rsaRes$Gene_ID[rsaRes$OPI_Rank < 999999])), collapse = ""));
  message(paste("Total #Hit Wells = ", sum(rsaRes$OPI_Rank < 999999), collapse = ""));
}

#' @keywords internal
handleOneGroup <- function(i, dataset, optsb, t.rank = NULL) {
  #- i is a vector with the position info in the sorted data.
  ## how many of them are up or low than cut-off.
  if (optsb$reverse) {
    i_max <- sum(dataset$Score[i] >= optsb$LB)
    i_min <- max(1, sum(dataset$Score[i] >= optsb$UB))
  } else {
    i_max <- sum(dataset$Score[i] <= optsb$UB)
    i_min <- max(1, sum(dataset$Score[i] <= optsb$LB))
  }

  ## t.r is the true rank, instead of order.
  if (is.null(t.rank)) {
    t.r <- i
  } else {
    t.r <- t.rank[i]
  }

  r <- OPIScore(t.r, nrow(dataset), i_min, i_max, optsb$bonferroni)

  return(cbind(LogP = r["logp"], OPI_Hit = as.numeric(seq(length(i)) <= r["cutoff"]), "#hitWell" = i_max, "#totalWell" = length(i), rank = i))
}

#' @keywords internal
OPIScore <- function(I_rank, N, i_min = 1, i_max = -1, bonferroni = FALSE) {
  #- Number of total siRNA for a gene.
  #- I_rank has the rank info of the siRNA to a gene.
  n_drawn <- length(I_rank)

  if (i_max == -1) {
    i_max <- n_drawn
  }

  r1 <- c(logp = 1, cutoff = 0)

  if (i_max < i_min) return(r1)

  best_logp <- 1
  cutoff    <- 0

  for (i in i_min:i_max) {
    #- Dealing with tie?
    if (i < i_max && I_rank[i] == I_rank[i + 1]) {
      next
    }

    #- i - 1, overlapping; I_rank[i], white ball; N - I_rank[i], black balls.
    logp <- phyper(i - 1, I_rank[i], N - I_rank[i], n_drawn, lower.tail = FALSE, log.p = TRUE)
    logp <- max(logp / log(10), -100)

    if (logp <= best_logp) {
      best_logp <- logp
      cutoff    <- i
    }
  }

  if (bonferroni) {
    best_logp <- best_logp + log(i_max - i_min + 1) / log(10)
  }

  return (c(logp = best_logp, cutoff = cutoff))
}

#' @keywords internal
OPI <- function(Groups, Scores, opts, Data = NULL) {
  t <- data.frame(cbind(Gene_ID = Groups, Score = Scores))
  Sorted_Order <- order(t$Score, decreasing = opts$reverse)
  Data <- Data[Sorted_Order, ]
  t <- t[Sorted_Order, ]

  ## get the ranks, "max" for the tie.
  t.rank <- rank(t$Score, ties.method = "max")
  #- Same genes are extracted in tapply.
  t <- do.call("rbind", tapply(seq(nrow(t)), list(t$Gene_ID), handleOneGroup, dataset = t, opts, t.rank))
  t <- cbind(Data, t[order(t[, "rank"]), ])

  # add OPI_Rank
  t <- t[order(t$LogP, t$Score * ifelse(opts$reverse, -1, 1)), ]
  t$OPI_Rank <- cumsum(t$OPI_Hit)
  t$OPI_Rank[t$OPI_Hit == 0] <- 999999

  # add Cutoff_Rank
  t <- t[order(t$Score * (ifelse(opts$reverse, -1, 1)), t$LogP), ]

  if (opts$reverse) {
    tmp <- (t$Score >= opts$LB)
  } else {
    tmp <- (t$Score <= opts$UB)
  }

  t$Cutoff_Rank <- cumsum(tmp)
  t$Cutoff_Rank[!tmp] <- 999999

# add EXP_Rank
  t$EXP_Rank <- pmin(t$OPI_Rank, t$Cutoff_Rank)
  t$EXP_Rank <- pmin(t$OPI_Rank, t$Cutoff_Rank)

  if (opts$reverse) {
    return(t[order(t$OPI_Rank, -t$Score), ])
  } else {
    return(t[order(t$OPI_Rank, t$Score), ])
  }
}

#' @keywords internal
.ff_rsaRatio <- function(masterPlate, dta, treatment, control, normMethod) {
  norm_res      <- .ff_masterPlateValue(masterPlate, dta, treatment, control, normMethod = normMethod)
  masterp_dta   <- norm_res[[1]]
  n_treat_plate <- norm_res[[2]]

  plate_ratio <- apply(masterp_dta, 1, .ff_ratio_rsa, n_treat_plate, n_treat_plate + 1)
  genes       <- sapply(as.character(names(plate_ratio)), function(x) unlist(strsplit(x," "))[[1]])
  rsaScore    <- data.frame(Gene_ID = genes, Well_ID = names(plate_ratio), Score = plate_ratio)

  rownames(rsaScore) <- NULL
  return(rsaScore)
}

#' @keywords internal
.ff_ratio_rsa <- function(x, y, z) {
  x <- as.numeric(x)
  a <- length(x)
  round(median(x[1:y], na.rm = TRUE) / median(x[z:a], na.rm = TRUE), 4)
}
