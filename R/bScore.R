#' Calculate B-score
#'
#' Calculate the B-score for plates belonging to the same master plate.
#' Positive / negative controls are removed from the calculation.
#'
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
#' bscore_res <- sapply(unique(exampleDat$MASTER_PLATE), bScore,
#'   exampleDat, treatment = "treatment", control = "control", simplify = FALSE)
#' @export
bScore <- function(masterPlate, dta, treatment, control, outFile = FALSE) {
  contPlate  <- dta[MASTER_PLATE == masterPlate & EXPERIMENT_MODIFICATION == control, unique(PLATE)]
  treatPlate <- dta[MASTER_PLATE == masterPlate & EXPERIMENT_MODIFICATION == treatment, unique(PLATE)]

  b_cont  <- sapply(contPlate, .ff_bscorePlate, dta)
  b_treat <- sapply(treatPlate, .ff_bscorePlate, dta)

  b_cominbed <- NA

  if (! identical(rownames(b_treat), rownames(b_cont))) {
    message("(!!) Differernt order between and treatment and control in master plate: ", masterPlate, "\n")
  } else {
    b_cominbed <- cbind(b_treat, b_cont)
  }

  if (outFile && !is.na(b_cominbed)) {
    write.table(b_cominbed, paste0(masterPlate, ".bscore.csv"), sep = "\t", quote = FALSE)
  }

  return(b_cominbed)
}

#' @keywords internal
.ff_bscorePlate <- function(plateName, dta) {
  ## internal function to calculate the B score;
  ## plateName, plate to be normalized;

  message("(II) Processing PLATE:", plateName, "\n")
  one_plate <- dta[PLATE == plateName & EXPERIMENT_TYPE == "sample"]

  norm_res <- dcast(one_plate, ROW_NAME ~ COL_NAME, value.var = "READOUT") %>%
    as.data.frame %>%
    set_rownames(as.character(.$ROW_NAME)) %>%
    inset("ROW_NAME", value = NULL) %>%
    as.matrix %>%
    medpolish(na.rm = TRUE, maxiter = 100)

  Bscore <- (norm_res$residuals / mad(norm_res$residuals, na.rm = TRUE)) %>% c

  well_name <- dcast(one_plate, ROW_NAME ~ COL_NAME, value.var = "WELL_CONTENT_NAME") %>%
    as.data.frame %>%
    set_rownames(as.character(.$ROW_NAME)) %>%
    inset("ROW_NAME", value = NULL) %>%
    as.matrix %>%
    c

  if (sum(is.na(well_name)) == 0) {
    names(Bscore) <- well_name
  }
  else {
    message("(!!) Wrong well names in plate: ", plateName)
  }

  return(Bscore)
}
