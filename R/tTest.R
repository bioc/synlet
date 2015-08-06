#' student's t-test basing on B-score
#'
#' Select hits by student's t-test using B-score from treatment and control plates.
#' @import magrittr
#' @import ggplot2
#' @param masterPlate the master plate to be analyzed
#' @param bScore normalized bScore
#' @param numTreat number of treatment plates
#' @param numCont number of control plates
#' @return A list containing student's t-test for each master plate
#' \itemize{
#'   \item pvalue: p-value of the t-test
#'   \item Treat_Cont: difference in bscore: treatment - control 
#'   \item p_adj: BH adjusted p-value
#' }
#' @references
#' Birmingham, A. et al. Statistical methods for analysis of high-throughput RNA interference screens. Nat Methods 6, 569-575 (2009).
#' @examples
#' bscore.res <- sapply(as.character(unique(exampleDat$MASTER_PLATE)), bScore,
#'   exampleDat, control = "control", treatment = "treatment", simplify = FALSE)
#' bscore.ttest  <- sapply(names(bscore.res), tTest, bscore.res, numTreat = 3,
#'   numCont = 3, simplify = FALSE, USE.NAMES = TRUE)
#' bscore.combined <- data.frame(do.call(rbind, lapply(names(bscore.ttest),
#'   function(x) if (!is.null(bscore.ttest[[x]])) {data.frame(MASTER_PLATE = x,
#'   siRNAs = rownames(bscore.ttest[[x]]), bscore.ttest[[x]])})))
#' @export
tTest <- function(masterPlate, bScore, numTreat, numCont) {
  # browser()
  message("Processing MASTER PLATE:", masterPlate,"\n")
  tem.da <- bScore[[masterPlate]]
  ## remove the "empty" wells.
  tem.da <- tem.da[grep("empty", rownames(tem.da), invert = TRUE), ]

  ttest.res <- t(apply(tem.da, 1, .ff_ttest, numTreat, numCont))
  ## NA rows are removed;
  tem.na <- ttest.res[is.na(ttest.res[, 1]), ]
  if (dim(tem.na)[1] != 0) {
    message("--- siRNAs with NA value in the t-test in master plate",
      masterPlate, "double check ---- ", "\n")
    for (siRNA in rownames(tem.na)) {
      message(siRNA, "\n")
    }
  }

  ttest.res <- ttest.res[!is.na(ttest.res[, 1]), ]

  if(dim(ttest.res)[1] != 0) {
    tem.1 <- p.adjust(ttest.res[, 1], "BH")
    ttest.res <- cbind(ttest.res, p_adj = tem.1)
    ttest.res
  } else {
    message("--- Warning: No valid t-test in masterplate", masterPlate, "---\n")
    }
}

#' @keywords internal
.ff_ttest <- function(x, numTreat, numCont) {
  # browser()
  state <- factor(rep(c("treatment", "control"),
    times = c(numTreat, numCont)),
    levels = c("treatment", "control"))
  x.da <- data.frame(value = as.numeric(x), state = state)
  ttest.res <-try(t.test(value ~ state, data = x.da), silent = TRUE)
  if (is(ttest.res,"try-error")) {
    tem.p <- NA
    tem.diff <- NA
  } else {
    tem.p <- ttest.res$"p.value"
    x.mean <- tapply(x.da$value, x.da$state, mean, na.rm = TRUE)
    tem.diff <- x.mean[1] - x.mean[2]
  }
  tem.1 <- c(tem.p, tem.diff)
  names(tem.1) <- c("pvalue", "Treat_Cont")
  tem.1
}
