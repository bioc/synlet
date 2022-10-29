#' @keywords internal
# .ff_plateNorm <- function(masterPlate, dat) {
#   # browser()
#   masterp.da <- dplyr::filter(dat, MASTER_PLATE == masterPlate,
#     EXPERIMENT_TYPE == "sample", WELL_CONTENT_NAME != "empty")
#   masterp.da$PLATE <- as.character(masterp.da$PLATE)
#   tem.1 <- doBy::summaryBy(READOUT ~ PLATE, data = masterp.da,
#     FUN = median, na.rm = TRUE)
#   colnames(tem.1)[2] <- "PLATE_MED"
#   masterp.da.1 <- merge(masterp.da, tem.1, by.x = "PLATE", by.y = "PLATE")
#   tem.2 <- masterp.da.1$READOUT / masterp.da.1$PLATE_MED
#   masterp.da.1 <- data.frame(masterp.da.1, NORMEDV = tem.2)
#   # masterp.da.1
# }

#' @keywords internal
.ff_plateNorm <- function(masterPlate, dta) {
  masterp_dta <- dta[MASTER_PLATE == masterPlate & EXPERIMENT_TYPE == "sample" & WELL_CONTENT_NAME != "empty"]

  masterp_norm <- masterp_dta[, .(PLATE_MED = median(READOUT, na.rm = TRUE)), by = "PLATE"] %>%
    merge(masterp_dta, ., by.x = "PLATE", by.y = "PLATE") %>%
    .[, NORMEDV := READOUT / PLATE_MED]

  return(masterp_norm)
}

#' @keywords internal
.ff_contsiRNANorm <- function(masterPlate, dta, contsiRNA) {
  masterp_dta <- dta[MASTER_PLATE == masterPlate & EXPERIMENT_TYPE == "sample" & WELL_CONTENT_NAME != "empty"]
  cont_dta    <- dta[MASTER_PLATE == masterPlate & WELL_CONTENT_NAME == contsiRNA]

  masterp_norm <- cont_dta[, .(CONTROL_MED = median(READOUT, na.rm = TRUE)), by = "PLATE"] %>%
    merge(masterp_dta, ., by.x = "PLATE", by.y = "PLATE", all.x = TRUE) %>%
    .[, NORMEDV := READOUT / CONTROL_MED]

  return(masterp_norm)
}

# #' @keywords internal
# .ff_contsiRNANorm <- function(masterPlate, dat, contsiRNA) {
#   # browser()
#   masterp.da <- dplyr::filter(dat, MASTER_PLATE == masterPlate,
#     EXPERIMENT_TYPE == "sample", WELL_CONTENT_NAME != "empty")
#   control.da <- dplyr::filter(dat, MASTER_PLATE == masterPlate,
#     WELL_CONTENT_NAME == contsiRNA)

#   masterp.da$PLATE <- as.character(masterp.da$PLATE)
#   control.da$PLATE <- as.character(control.da$PLATE)
#   tem.1 <- doBy::summaryBy(READOUT ~ PLATE, data = control.da, FUN = median, na.rm = TRUE)
#   # if (! all(complete.cases(tem.1))) browser()
#   colnames(tem.1)[2] <- "CONTROL_MED"
#   ## contain all masterplate data, even in the condition no control available.
#   masterp.da.1 <- merge(masterp.da, tem.1, by.x = "PLATE", by.y = "PLATE", all.x = TRUE)
#   ## normalize to the median.
#   tem.2 <- masterp.da.1$READOUT / masterp.da.1$CONTROL_MED
#   masterp.da.1 <- cbind(masterp.da.1, NORMEDV = tem.2)
#   # masterp.da.1
# }

#' @keywords internal
.ff_masterPlateValue <- function(masterPlate, dat, treatment, control, normMethod = "PLATE") {
  if (normMethod == "PLATE") {
    masterp.da.1 <- .ff_plateNorm(masterPlate, dat)
  } else {
    masterp.da.1 <- .ff_contsiRNANorm(masterPlate, dat, normMethod)
  }

  treat.p <- dplyr::filter(masterp.da.1, EXPERIMENT_MODIFICATION == treatment) %>%
    extract2("PLATE") %>%
    as.character %>%
    unique

  cont.p <- dplyr::filter(masterp.da.1, EXPERIMENT_MODIFICATION == control) %>%
    extract2("PLATE") %>%
    as.character %>%
    unique

  # browser()
  if (length(treat.p) != 0 & length(cont.p) != 0) {
    ## treatment first, followed by controls;
    plate.paired <- c(treat.p, cont.p)
    masterp.da.2 <- reshape2::dcast(masterp.da.1, WELL_CONTENT_NAME ~ PLATE,
      value.var = "NORMEDV")
    rownames(masterp.da.2) <- masterp.da.2$WELL_CONTENT_NAME
    masterp.da.2 <- as.data.frame(masterp.da.2[, -1])
    masterp.da.2 <- masterp.da.2[, plate.paired]
    master.norm <- list(normedV = masterp.da.2, treat.num = length(treat.p),
      cont.num = length(cont.p))
    # master.norm

  } else {
    stop("--- Empty control or treament plates in master plate:",
      masterPlate, "---\n")
  }
}
