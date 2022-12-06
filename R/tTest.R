#' student's t-test on B-score
#'
#' Select hits by student's t-test using B-score from treatment and control plates.
#'
#' @param mtx b-score matrix.
#' @param n_treat number of treatment plates
#' @param n_cont number of control plates
#' @return A list containing student's t-test for each master plate
#' \itemize{
#'   \item pvalue: p-value of the t-test
#'   \item Treat_Cont: difference in bscore: treatment - control
#'   \item p_adj: BH adjusted p-value
#' }
#' @references
#' Birmingham, A. et al. Statistical methods for analysis of high-throughput RNA interference screens. Nat Methods 6, 569-575 (2009).
#' @examples
#' data(example_dt)
#' bscore_res <- sapply(unique(example_dt$MASTER_PLATE), bScore,
#'   example_dt, control = "control", treatment = "treatment", simplify = FALSE)
#' tTest(bscore_res$P001, 3, 3)
#' @export
tTest <- function(mtx, n_treat, n_cont) {
  mtx <- mtx[grep("empty", rownames(mtx), invert = TRUE), ]

  res    <- t(apply(mtx, 1, .ff_ttest, n_treat, n_cont))
  res_na <- res[is.na(res[, "pvalue"]), ]

  if (nrow(res_na) != 0) {
    message("(II) NA generated for following siRNAs in t-test.")
    for (siRNA in rownames(res_na)) {
      message(siRNA, "\n")
    }
  }

  res <- res[!is.na(res[, "pvalue"]), ]

  #- replace number with column names.
  if (nrow(res) != 0) {
    res <- cbind(res, p_adj = p.adjust(res[, "pvalue"], "BH"))
  } else {
    message("(II) No valid t-test.")
  }

  return(res)
}

#' @keywords internal
.ff_ttest <- function(x, n_treat, n_cont) {
  # browser()
  cond <- factor(rep(c("treatment", "control"), times = c(n_treat, n_cont)), levels = c("treatment", "control"))
  dta  <- data.frame(value = x, cond = cond)

  res <- try(t.test(value ~ cond, data = dta), silent = TRUE)

  if (is(res, "try-error")) {
    p <- NA
    v <- NA
  } else {
    p <- res$"p.value"
    m <- tapply(dta$value, dta$cond, mean, na.rm = TRUE)
    v <- m[1] - m[2]
  }

  t_res <- c(p, v) %>% set_names(c("pvalue", "Treat_Cont"))
  return(t_res)
}
