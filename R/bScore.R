#' Calculate B-score
#'
#' Calculate the B-score for plates belonging to the same master plate.
#' Positive / negative controls are removed from the calculation.
#'
#' @import magrittr
#' @import ggplot2
#' @importFrom dplyr filter
#' @importFrom reshape2 dcast
#' @param masterPlate a maste plate to be normalized.
#' @param dta synthetic lethal RNAi screen data.
#' @param treatment the treatment experiment condition in EXPERIMENT_MODIFICATION
#' @param control the control experiment condition in EXPERIMENT_MODIFICATION.
#' @param outFile should calculated B-score files be written to the current folder? File names is (masterPlate).bscore.csv.
#' @return A list contains B-score for each master plate, treatment plates are the
#'   first columns, followed by control plates
#' @references
#' Brideau, C., Gunter, B., Pikounis, B. & Liaw, A. Improved statistical methods for hit selection in high-throughput screening. J. Biomol. Screen. 8, 634-647 (2003).
#' @examples
#' bscore.res <- sapply(as.character(unique(exampleDat$MASTER_PLATE)), bScore,
#'   exampleDat, control = "control", treatment = "treatment", simplify = FALSE)
#' @export
bScore <- function(masterPlate, dta, treatment, control, outFile = FALSE) {
  #- Extract contral plates
  # contPlate <- dplyr::filter(dta, EXPERIMENT_MODIFICATION == control,
  #   MASTER_PLATE == masterPlate) %>%
  # extract(, 1) %>%
  # unique %>%
  # as.character

  #- Tested.
  contPlate <- dta[MASTER_PLATE == masterPlate & EXPERIMENT_MODIFICATION == control, unique(PLATE)]

  # treatPlate <- dplyr::filter(dta, EXPERIMENT_MODIFICATION == treatment,
  #   MASTER_PLATE == masterPlate) %>%
  # extract(, 1) %>%
  # unique %>%
  # as.character

  #- Tested.
  treatPlate <- dta[MASTER_PLATE == masterPlate & EXPERIMENT_MODIFICATION == treatment, unique(PLATE)]

  b_cont  <- sapply(contPlate, .ff_bscorePlate, dta)
  b_treat <- sapply(treatPlate, .ff_bscorePlate, dta)

  b_cominbed <- NA

  if (! identical(rownames(b_treat), rownames(b_cont))) {
    message("(!!) Differernt order between and treatment and control in master plate: ", masterPlate, "\n")
    break()
  } else {
    b_cominbed <- cbind(b_treat, b_cont)
  }

  # browser()
  if (outFile && !is.na(b_cominbed)) {
    write.table(b_cominbed, paste0(masterPlate, ".bscore.csv"), sep = "\t", quote = FALSE)
  }

  return(b_cominbed)
}

#' @keywords internal
.ff_bscorePlate <- function(plateName, dta) {
  ## internal function to calculate the B score;
  ## plateName, plate to be normalized;
  ## dta, screen data;
  # browser()

  #- !! dta is a data.table.
  message("(II) Processing PLATE:", plateName, "\n")
  # plateName <- as.numeric(plateName)
  # tem.1 <- dplyr::filter(dta, PLATE == as.numeric(plateName), EXPERIMENT_TYPE == "sample")
  one_plate <- dta[PLATE == plateName & EXPERIMENT_TYPE == "sample"]

  # master_name <- tem.1[1, "MASTER_PLATE"] %>% as.character
  # master_name <- one_plate$MASTER_PLATE[1]

  # tem.3 <- reshape2::dcast(tem.1, ROW_NAME ~ COL_NAME, value.var = "READOUT")
  # rownames(tem.3) <- tem.3[, 1]
  # tem.3 <- as.matrix(tem.3[, -1])

  #- !! Tested.
  # cell_readout <- dcast(one_plate, ROW_NAME ~ COL_NAME, value.var = "READOUT") %>%
  #   as.data.frame %>%
  #   set_rownames(as.character(.$ROW_NAME)) %>%
  #   inset("ROW_NAME", value = NULL) %>%
  #   as.matrix

  norm_res <- dcast(one_plate, ROW_NAME ~ COL_NAME, value.var = "READOUT") %>%
    as.data.frame %>%
    set_rownames(as.character(.$ROW_NAME)) %>%
    inset("ROW_NAME", value = NULL) %>%
    as.matrix %>%
    medpolish(na.rm = TRUE, maxiter = 100)

  ## get the residuals.
  # tem.4 <- medpolish(tem.3, na.rm = TRUE, maxiter = 100)
  # norm_res <- medpolish(cell_readout, na.rm = TRUE, maxiter = 100)

  ## bscore, mad in R is a little different from the original paper with a coefficent.
  #- !! Here.
  # Bscore <- tem.4$residuals / mad(tem.4$residuals, na.rm = TRUE)
  # Bscore <- c(Bscore)

  Bscore <- (norm_res$residuals / mad(norm_res$residuals, na.rm = TRUE)) %>% c
  # Bscore <- c(Bscore)

  # tem.6 <- reshape2::dcast(tem.1, ROW_NAME ~ COL_NAME,
  #   value.var = "WELL_CONTENT_NAME")
  # rownames(tem.6) <- tem.6[, 1]
  # tem.6 <- tem.6[, -1]
  # tem.6 <- as.character(unlist(tem.6))

  well_name <- dcast(one_plate, ROW_NAME ~ COL_NAME, value.var = "WELL_CONTENT_NAME") %>%
    as.data.frame %>%
    set_rownames(as.character(.$ROW_NAME)) %>%
    inset("ROW_NAME", value = NULL) %>%
    as.matrix %>%
    c

  # all.equal(tem.6, test.6)

  if (sum(is.na(well_name)) == 0) {
    # names(Bscore) <- tem.6
    names(Bscore) <- well_name
  }
  else {
    message("(!!) Wrong well names in plate: ", plateName)
  }

  # Bscore <- sort(Bscore)
  # Bscore_2 <- sort(Bscore_2)
  # all.equal(Bscore_2, Bscore)

  return(Bscore)
}
