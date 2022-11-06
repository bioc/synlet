#' Select hits basing on median +- k*MAD
#'
#' Select hits basing on median +- k*MAD, by default k is three.
#' @param masterPlate the master plate to analysis
#' @param dat synthetic lethal RNAi screen data
#' @param k cutoff for selecting hits, default is three
#' @param treatment the treatment condition in EXPERIMENT_MODIFICATION
#' @param control the control condition in EXPERIMENT_MODIFICATION
#' @param outFile whether or not write the median normalized results
#' @param normMethod normalization methods to be used. If "PLATE", the raw readouts are normalized by plate median, otherwise use median provided control siRNA.
#' @return A data.frame contains the hits selection results.
#' \itemize{
#'   \item MASTER_PLATE: location of siRNA
#'   \item treat_cont_ratio: ratio of treatment / control
#'   \item treat_median: median value of treatment plates
#'   \item control_median: median value of control plates
#'   \item Hits: Is this siRNA a hit?
#' }
#' @references
#' Chung,N.etal. Median absolute deviation to improve hits election for genome-scale RNAi screens. J. Biomol. Screen. 13, 149-158 (2008).
#' @examples
#' madSelection <- sapply(as.character(unique(exampleDat$MASTER_PLATE)),
#'   madSelect, exampleDat, control = "control",
#'   treatment = "treatment", simplify = FALSE)
#' madSelection.c <- do.call(rbind,
#'   lapply(names(madSelection), function(x) madSelection[[x]]))
#' @export
madSelect <- function(masterPlate, dat, k = 3, treatment, control, outFile = FALSE, normMethod = "PLATE") {
  norm_res      <- .ff_masterPlateValue(masterPlate, dat, treatment, control, normMethod = normMethod)
  masterp_dta   <- norm_res[[1]]
  n_treat_plate <- norm_res[[2]]

  #- All plate belongs to the same master plate are concatenate together.
  #- !! I would like drop append == TRUE as only one master plate exist
  if (outFile) {
    write.table(masterp_dta, paste0(masterPlate, ".median.norm.txt"), quote = FALSE, sep = "\t")
  }

  #- Get confused here.
  plate_ratio    <- t(apply(masterp_dta, 1, .ff_ratio_madS, n_treat_plate, n_treat_plate + 1))
  plate_r_med    <- median(plate_ratio[, 1])
  Hits           <- rep("No", nrow(plate_ratio))
  Hit_sele       <- (plate_ratio[, 1] < (plate_r_med - k * mad(plate_ratio[, 1])))
  Hits[Hit_sele] <- "Yes"

  res <- data.frame(MASTER_PLATE = masterPlate, plate_ratio, Hits = Hits)
  # masterp.da.ratio
}

#' @keywords internal
.ff_ratio_madS <- function(x, y, z) {
  a           <- length(x)
  med_treat   <- median(x[1:y], na.rm = TRUE)
  med_cont    <- median(x[z:a], na.rm = TRUE)
  tr_co_ratio <- med_treat / med_cont

  l <- c(tr_co_ratio, med_treat, med_cont)
  names(l) <- c("treat_cont_ratio", "treat_median", "control_median")
  return(l)
}
