#' Select hits by the rank product method
#'
#' Select hits by rank product methods by comparing treatment and control.
#' @import magrittr
#' @import ggplot2
#' @importFrom RankProd RP
#' @param masterPlate the master plate to be analyzed
#' @param dat synthetic lethal RNAi screen data
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
#'
#' Hong, F. et al. RankProd: a bioconductor package for detecting differentially expressed genes in meta-analysis. Bioinformatics 22, 2825-2827 (2006).
#' @examples
#' rankp.res <- sapply(as.character(unique(exampleDat$MASTER_PLATE)),
#'   rankProdHits, exampleDat, control = "control", treatment = "treatment",
#'   simplify = FALSE)
#' rankp.c <- data.frame(do.call(rbind,
#'   lapply(names(rankp.res), function(x) rankp.res[[x]])))
#' @export
rankProdHits <- function(masterPlate, dat, treatment, control,
  normMethod = "PLATE") {
  # browser()
  if(class(dat) != "data.frame") stop("incorrect input!")
  master.norm <- .ff_masterPlateValue(masterPlate, dat, treatment, control,
    normMethod = normMethod)
  masterp.da <- master.norm[[1]]
  treat.num <- master.norm[[2]]
  cont.num <- master.norm[[3]]

  siRNA.cl <- rep(c(0, 1), times = c(treat.num, cont.num))
  siRNA.rp <- RankProd::RP(masterp.da, siRNA.cl, logged = FALSE)
  masterp.rp <- data.frame(MASTER_PLATE = masterPlate,
    pvalue_treat_lowerthan_cont = siRNA.rp$pval[, 1],
    FDR_treat_lowerthan_cont = siRNA.rp$pfp[, 1],
    treat_cont_log2FC = log2(siRNA.rp$AveFC[, 1]))
  # masterp.rp
}
