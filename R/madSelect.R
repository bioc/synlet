#' Select hits basing on median +- k*MAD
#'
#' Select hits basing on median +- k*MAD, by default k is three.
#' @import magrittr
#' @import ggplot2
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
#' Chung,N.etal. Medianabsolutedeviationtoimprovehitselectionforgenome- scale RNAi screens. J. Biomol. Screen. 13, 149-158 (2008).
#' @examples
#' madSelection <- sapply(as.character(unique(exampleDat$MASTER_PLATE)),
#'   madSelect, exampleDat, control = "control",
#'   treatment = "treatment", simplify = FALSE)
#' madSelection.c <- do.call(rbind,
#'   lapply(names(madSelection), function(x) madSelection[[x]]))
#' @export
madSelect <- function(masterPlate, dat, k = 3, treatment, control,
   outFile = FALSE, normMethod = "PLATE") {
  cutoff = k
  if(class(dat) != "data.frame") stop("incorrect input!")
  master.norm <- .ff_masterPlateValue(masterPlate, dat, treatment, control,
    normMethod = normMethod)
  masterp.da <- master.norm[[1]]
  treat.num <- master.norm[[2]]
  cont.num <- master.norm[[3]]

  if (outFile) {
    write.table(masterp.da, paste0(masterPlate, ".median.norm.da.txt"),
      append = TRUE, quote = FALSE, sep = "\t",
      row.names = FALSE, col.names = FALSE)
  }

  plate.ratio <- t(apply(masterp.da, 1, .ff_ratio_madS, treat.num, treat.num + 1))
  plate.r.median <- median(plate.ratio[, 1])
  # browser()
  Hits <- rep("No", length(plate.ratio[, 1]))
  Hit.sele <- (plate.ratio[,1] < (plate.r.median - cutoff * mad(plate.ratio[,1])))
  Hits[Hit.sele] <- "Yes"

   masterp.da.ratio <- data.frame(MASTER_PLATE = masterPlate, plate.ratio, Hits = Hits)
  # masterp.da.ratio
}


#' @keywords internal
.ff_ratio_madS <- function(x, y, z) {
  x <- as.numeric(x)
  a <- length(x)
  tem.1 <- round(median(x[1:y], na.rm = TRUE) / median(x[z:a], na.rm = TRUE), 4)
  tem.treat <- median(x[1:y], na.rm = TRUE)
  tem.cont <- median(x[z:a], na.rm = TRUE)

  tem.2 <- c(tem.1, tem.treat, tem.cont)
  names(tem.2) <- c("treat_cont_ratio", "treat_median", "control_median")
  tem.2
}
