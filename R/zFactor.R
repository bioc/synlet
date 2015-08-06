#' Calcualte the Z and Z' factor
#'
#' calcualte the Z and Z' factor for each plate.
#' @import magrittr
#' @import ggplot2
#' @importFrom dplyr filter
#' @importFrom doBy summaryBy
#' @importFrom reshape2 dcast
#' @param dat synthetic lethal RNAi screen data.
#' @param negativeCon the negative control used in the WELL_CONTENT_NAME.
#' @param positiveCon the positive control used in the WELL_CONTENT_NAME.
#' @param useMean use mean to calcualate z factor and z' factor by default; otherwise use median.
#' @examples
#' zFactor(exampleDat, negativeCon = "scrambled control si1", positiveCon = "PLK1 si1")
#' @return  A data.frame contains z factor and z' factor
#' @references
#' Zhang J.H., Chung T.D. & Oldenburg K.R. A simple statistical parameter for use in evaluation and validation of high throughput screening assays. J. Biomol. Screen. B, 4 67-73 (1999).
#' Birmingham,A. et al. (2009) Statistical methods for analysis of high-throughput RNA interference screens. Nat Methods, 6, 569-575.
#' @export
zFactor <- function(dat, negativeCon, positiveCon, useMean = TRUE) {
  if(class(dat) != "data.frame") stop("incorrect input!")
  tem.1 <- as.character(dat[, "WELL_CONTENT_NAME"])
  tem.1[dat$EXPERIMENT_TYPE == "sample"] <- "treat"
  # browser()
  tem.1[grep(negativeCon, tem.1)] <- "neg_control"
  tem.1[grep(positiveCon, tem.1)] <- "posi_control"
  data.f3 <- cbind(dat, condition = tem.1)
  data.f3$PLATE <- factor(data.f3$PLATE)
  data.f3 <- dplyr::filter(data.f3,
    condition == "treat" | condition == "neg_control" | condition == "posi_control")
  data.f3$condition <- factor(data.f3$condition)
  if (useMean) {
    data.f4 <- doBy::summaryBy(READOUT ~ PLATE + condition + EXPERIMENT_MODIFICATION,
      data = data.f3, FUN = c(sd, mean), na.rm = TRUE)
  }
  else {
    data.f4 <- doBy::summaryBy(READOUT ~ PLATE + condition + EXPERIMENT_MODIFICATION,
      data = data.f3, FUN = c(sd, median), na.rm = TRUE)
  }

  colnames(data.f4)[4:5] <- c("sd_sample","m_sample")
  data.f5 <- reshape2::dcast(data.f4, PLATE ~ condition, value.var = "sd_sample")
  colnames(data.f5)[2:4] <- paste("sd", colnames(data.f5)[2:4],sep = "_")
  data.f6 <- reshape2::dcast(data.f4, PLATE ~ condition, value.var = "m_sample")
  colnames(data.f6)[2:4] <- paste("m", colnames(data.f6)[2:4], sep = "_")

  if(identical(data.f5[, 1], data.f6[, 1])) {
    data.f7 <- cbind(data.f5, data.f6[, -1])
    rownames(data.f7) <- data.f7[, 1]
    data.f7 <- data.f7[complete.cases(data.f7), ]
    message("--- number of plates to calculate Z factor: ", dim(data.f7)[1],"---\n")
    tem.1 <- t(apply(data.f7, 1, .ff_z))
    colnames(tem.1) <- c("zFactor", "zPrimeFactor")
    tem.1
  } else {
    message("ERROR: mean/median and sd are not paired!","\n")
  }
}

#' @keywords internal
.ff_z <- function(x) {
  x <- as.numeric(x)
  sd_nega <- x[2]
  sd_posi <- x[3]
  sd_sample <- x[4]
  m_nega <- x[5]
  m_posi <- x[6]
  m_sam <- x[7]
  z_factor <- 1 - 3 * (sd_nega + sd_sample) / abs(m_nega - m_sam)
  z_primef <- 1 - 3 * (sd_nega + sd_posi) / abs(m_nega - m_posi)

  return(c(z_factor, z_primef))
}
