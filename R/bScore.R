#' Calculate B-score
#'
#' Calculate the B-score for plates belonging to the same master plate.
#' Positive / negative controls are removed from the calculation.
#'
#' @import magrittr
#' @import ggplot2
#' @importFrom dplyr filter
#' @importFrom reshape2 dcast
#' @param masterPlate a maste plate to be normalized
#' @param dat synthetic lethal RNAi screen data
#' @param treatment the treatment experiment condition in EXPERIMENT_MODIFICATION
#' @param control the control experiment condition in EXPERIMENT_MODIFICATION
#' @param outFile should calculated B-score files be written to the current folder? File names is (masterPlate).bscore.csv.
#' @return A list contains B-score for each master plate, treatment plates are the
#' first columns, followed by control plates
#' @references
#' Brideau, C., Gunter, B., Pikounis, B. & Liaw, A. Improved statistical methods for hit selection in high-throughput screening. J. Biomol. Screen. 8, 634-647 (2003).
#' @examples
#' bscore.res <- sapply(as.character(unique(exampleDat$MASTER_PLATE)), bScore,
#'   exampleDat, control = "control", treatment = "treatment", simplify = FALSE)
#' @export
bScore <- function(masterPlate, dat, treatment, control, outFile = FALSE) {
  # browser()
  if(class(dat) != "data.frame") stop("incorrect input!")
  contPlate <- dplyr::filter(dat, EXPERIMENT_MODIFICATION == control,
    MASTER_PLATE == masterPlate) %>%
  extract(, 1) %>%
  unique %>%
  as.character

  treatPlate <- dplyr::filter(dat, EXPERIMENT_MODIFICATION == treatment,
    MASTER_PLATE == masterPlate) %>%
  extract(, 1) %>%
  unique %>%
  as.character

  tem.treat <- sapply(treatPlate, .ff_bscorePlate, dat)
  tem.cont <- sapply(contPlate, .ff_bscorePlate, dat)

  if (! identical(rownames(tem.treat), rownames(tem.cont))) {
    message("Error, differernt order between and treatment and control", "\n")
  } else {
    tem.out <- cbind(tem.treat, tem.cont)
  }
  # browser()
  if (outFile) {
    write.table(tem.out, paste0(masterPlate, ".bscore.csv"), sep = "\t",
      quote = FALSE)
  }
  tem.out
}

#' @keywords internal
.ff_bscorePlate <- function(plateName, dat) {
  ## internal function to calculate the B score;
  ## plateName, plate to be normalized;
  ## dat, screen data;
  # browser()
  message("---Processing PLATE:", plateName, "---\n")
  plateName <- as.numeric(plateName)
  tem.1 <- dplyr::filter(dat, PLATE == plateName, EXPERIMENT_TYPE == "sample")
  master.n <- tem.1[1, "MASTER_PLATE"] %>% as.character

  tem.3 <- reshape2::dcast(tem.1, ROW_NAME ~ COL_NAME, value.var = "READOUT")
  rownames(tem.3) <- tem.3[, 1]
  tem.3 <- tem.3[, -1]
  ## get the residuals.
  tem.4 <- medpolish(tem.3, na.rm = TRUE, maxiter = 100)
  ## bscore, mad in R is a little different from the original paper with a coefficent.
  Bscore <- tem.4$residuals / mad(tem.4$residuals, na.rm = TRUE)
  Bscore <- c(Bscore)

  tem.6 <- reshape2::dcast(tem.1, ROW_NAME ~ COL_NAME,
    value.var = "WELL_CONTENT_NAME")
  rownames(tem.6) <- tem.6[, 1]
  tem.6 <- tem.6[, -1]
  tem.6 <- as.character(unlist(tem.6))

  if (sum(is.na(tem.6)) == 0) {
    names(Bscore) <- tem.6
  }
  else {
    message("Error: wrong well names in plate:", plateName)
  }
  Bscore
}
