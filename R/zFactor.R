#' Calcualte the Z and Z' factor
#'
#' calcualte the Z and Z' factor for each plate.
#'
#' @param dta synthetic lethal RNAi screen data.
#' @param negativeCon the negative control used in the WELL_CONTENT_NAME.
#' @param positiveCon the positive control used in the WELL_CONTENT_NAME.
#' @param useMean use mean to calcualate z factor and z' factor by default; otherwise use median.
#' @return  A data.frame contains z factor and z' factor
#' @examples
#' data(example_dt)
#' res <- zFactor(example_dt, negativeCon = "scrambled control si1", positiveCon = "PLK1 si1")
#' @references
#' Zhang J.H., Chung T.D. & Oldenburg K.R. A simple statistical parameter for use in evaluation and validation of high throughput screening assays. J. Biomol. Screen. B, 4 67-73 (1999).
#' Birmingham,A. et al. (2009) Statistical methods for analysis of high-throughput RNA interference screens. Nat Methods, 6, 569-575.
#' @export
zFactor <- function(dta, negativeCon, positiveCon, useMean = TRUE) {
  well_type <- WELL_CONTENT_NAME <- EXPERIMENT_TYPE <- READOUT <- PLATE <- NULL

  dta_2 <- copy(dta) %>%
    .[, well_type := WELL_CONTENT_NAME] %>%
    .[EXPERIMENT_TYPE == "sample", well_type := "treat"] %>%
    .[well_type %like% negativeCon, well_type := "neg_control"] %>%
    .[well_type %like% positiveCon, well_type := "posi_control"] %>%
    .[well_type %in% c("treat", "neg_control", "posi_control")]

  if (useMean) {
    dta_3 <- dta_2[, .(sd_value = sd(READOUT, na.rm = TRUE), m_value = mean(READOUT, na.rm = TRUE)), by = c("PLATE", "well_type", "EXPERIMENT_MODIFICATION")]
  } else {
    dta_3 <- dta_2[, .(sd_value = sd(READOUT, na.rm = TRUE), m_value = median(READOUT, na.rm = TRUE)), by = c("PLATE", "well_type", "EXPERIMENT_MODIFICATION")]
  }

  dta_sd <- dcast(dta_3, PLATE ~ well_type, value.var = "sd_value") %>% setnames(c(2:4), paste0("sd_", colnames(.)[2:4]))
  dta_m  <- dcast(dta_3, PLATE ~ well_type, value.var = "m_value") %>% setnames(c(2:4), paste0("m_", colnames(.)[2:4]))

  if(identical(dta_sd$PLATE, dta_m$PLATE)) {
    dta_all <- cbind(dta_sd, dta_m[, PLATE := NULL])
    mtx     <- as.matrix(dta_all[, PLATE := NULL]) %>% set_rownames(dta_sd$PLATE) %>% na.omit
    mtx     <- mtx[, c("sd_neg_control", "sd_posi_control", "sd_treat", "m_neg_control", "m_posi_control", "m_treat")]

    message("(II) Number of plates to calculate Z/Z' factor: ", nrow(mtx), "\n")

    res <- t(apply(mtx, 1, .ff_z)) %>% set_colnames(c("zFactor", "zPrimeFactor"))
  } else {
    message("(EE) Plate names for mean/median and sd are not paired!", "\n")
    res <- NULL
  }

  return(res)
}

#' @keywords internal
.ff_z <- function(x) {
  sd_nega  <- x[1]
  sd_posi  <- x[2]
  sd_treat <- x[3]
  m_nega   <- x[4]
  m_posi   <- x[5]
  m_treat  <- x[6]

  z_factor <- 1 - 3 * (sd_nega + sd_treat) / abs(m_nega - m_treat)
  z_primef <- 1 - 3 * (sd_nega + sd_posi) / abs(m_nega - m_posi)

  return(c(z_factor, z_primef))
}
