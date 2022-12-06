#' Select hits by the rank product method
#'
#' Select hits by rank product methods by comparing treatment and control.
#'
#' @param masterPlate the master plate to be analyzed
#' @param dta synthetic lethal RNAi screen data
#' @param treatment the treatment condition in EXPERIMENT_MODIFICATION
#' @param control the control condition in EXPERIMENT_MODIFICATION
#' @param normMethod normalization methods to be used. If "PLATE", the raw readouts are normalized by plate median, otherwise use provided control siRNA
#' @return A list contains results by the rank product method for each master plate.
#' \itemize{
#'   \item MASTER_PLATE: location of siRNA
#'   \item pvalue_treat_lowerthan_cont: p-value for the hypothesis that treatment has lower normalized readout compared to control
#'   \item FDR_treat_lowerthan_cont: FDR value
#'   \item treat_cont_log2FC: log2 fold change of treatment / control
#' }
#' @references
#' Breitling, R., Armengaud, P., Amtmann, A. & Herzyk, P. Rank products: a simple, yet powerful, new method to detect differentially regulated genes in replicated microarray experiments. FEBS Lett 573, 83-92 (2004).
#' Hong, F. et al. RankProd: a bioconductor package for detecting differentially expressed genes in meta-analysis. Bioinformatics 22, 2825-2827 (2006).
#' @examples
#' data(example_dt)
#' res <- sapply(unique(example_dt$MASTER_PLATE),
#'               rankProdHits,
#'               example_dt,
#'               control   = "control",
#'               treatment = "treatment",
#'               simplify  = FALSE)
#' @export
rankProdHits <- function(masterPlate, dta, treatment, control, normMethod = "PLATE") {
  norm_res      <- .ff_masterPlateValue(masterPlate, dta, treatment, control, normMethod = normMethod)
  masterp_dta   <- norm_res[[1]]
  n_treat_plate <- norm_res[[2]]
  n_cont_plate  <- norm_res[[3]]

  siRNA_cl   <- rep(c(0, 1), times = c(n_treat_plate, n_cont_plate))
  siRNA_rp   <- RankProd::RP(masterp_dta, siRNA_cl, logged = FALSE)
  masterp_rp <- data.frame(MASTER_PLATE                = masterPlate,
                           pvalue_treat_lowerthan_cont = siRNA_rp$pval[, 1],
                           FDR_treat_lowerthan_cont    = siRNA_rp$pfp[, 1],
                           treat_cont_log2FC           = log2(siRNA_rp$AveFC[, 1]))
  return(masterp_rp)
}
