#- norm_method is either "PLATE" or siRNA to be used (e.g., "lipid only" in WELL_CONTENT_NAME)
#' @keywords internal
.ff_norm <- function(masterPlate, dta, norm_method = "PLATE") {
  WELL_CONTENT_NAME <- MASTER_PLATE <- EXPERIMENT_TYPE <- READOUT <- NORMEDV <- BACKGROUND_MED <- NULL
  masterp_dta <- dta[MASTER_PLATE == masterPlate & EXPERIMENT_TYPE == "sample" & WELL_CONTENT_NAME != "empty"]

  if (norm_method == "PLATE") {
    norm_dta <- masterp_dta[, .(BACKGROUND_MED = median(READOUT, na.rm = TRUE)), by = "PLATE"]
  } else {
    norm_dta <- dta[MASTER_PLATE == masterPlate & WELL_CONTENT_NAME == norm_method] %>%
      .[, .(BACKGROUND_MED = median(READOUT, na.rm = TRUE)), by = "PLATE"]
  }

  masterp_norm <- merge(masterp_dta, norm_dta, by.x = "PLATE", by.y = "PLATE", all.x = TRUE) %>%
    .[, NORMEDV := READOUT / BACKGROUND_MED]

  return(masterp_norm)
}

#- normMethod is either "PLATE", or contron siRNA names. Former is prefered.
#' @keywords internal
.ff_masterPlateValue <- function(masterPlate, dta, treatment, control, normMethod = "PLATE") {
  EXPERIMENT_MODIFICATION <- PLATE <- NULL
  masterp_dta <- .ff_norm(masterPlate, dta, normMethod)

  treat_p <- masterp_dta[EXPERIMENT_MODIFICATION == treatment, unique(PLATE)]
  cont_p  <- masterp_dta[EXPERIMENT_MODIFICATION == control, unique(PLATE)]

  if (length(treat_p) != 0 & length(cont_p) != 0) {
    ## treatment first, followed by controls;
    plate_paired <- c(treat_p, cont_p)
    masterp_dta_w <- dcast(masterp_dta, WELL_CONTENT_NAME ~ PLATE, value.var = "NORMEDV") %>%
      as.data.frame %>%
      set_rownames(.$WELL_CONTENT_NAME) %>%
      inset("WELL_CONTENT_NAME", value = NULL) %>%
      extract(, plate_paired)

    master_norm <- list(normedV = masterp_dta_w, n_treat_plate = length(treat_p), n_cont_plate = length(cont_p))
  } else {
    stop("(EE) Empty control or treament plates in master plate:", masterPlate, "\n")
  }
}
