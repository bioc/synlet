#' Scatter plot of RNAi screen results
#'
#' Produce a single plot for readous of each plate, with the option of highlighting
#' specific signals, like positive/negative controls.
#'
#' @import magrittr
#' @import ggplot2
#' @param dta synthetic lethal RNAi screen data
#' @param scatter_colour colour for different signals
#' @param controlOnly whether or not to plot control wells only
#' @param control_name names of control siRNAs.
#' @return a ggplot object
#' @examples
#' scatterPlot(exampleDat, controlOnly = FALSE, colour = rainbow(10),
#'   "PLK1 si1", "scrambled control si1", "lipid only")
#' @export
# scatterPlot <- function(dta, colour = rainbow(10), controlOnly = FALSE, ...) {
scatterPlot <- function(dta, scatter_colour = rainbow(10), controlOnly = FALSE, control_name = NULL) {
  dta_2 <- copy(dta)

  if (controlOnly == FALSE) {
    if (is.null(control_name)) {
      message("(II) Better to provide names for positive/negative controls.")
      dta_2[, new_cond := "knock-down signals"]
    } else {
      dta_2[, new_cond := WELL_CONTENT_NAME] %>%
        .[!(WELL_CONTENT_NAME %in% control_name), new_cond := "knock-down signals"]
    }
  } else {
    if (is.null(control_name)) {
      stop("(EE) Please provide names for positive/negative controls if controlOnly is TRUE.")
    } else {
      dta_2 <- dta_2[WELL_CONTENT_NAME %in% control_name][, new_cond := WELL_CONTENT_NAME]
    }
  }

  dta_2 <- dta_2[order(EXPERIMENT_MODIFICATION)] %>%
    .[, new_name := paste(PLATE, EXPERIMENT_MODIFICATION, sep = " ") %>% factor]

  ggplot(as.data.frame(dta_2), aes(new_name, READOUT, colour = new_cond)) +
    geom_point() +
    theme_bw(14) +
    theme(panel.grid.major     = element_blank(),
          panel.grid.minor     = element_blank(),
          panel.background     = element_rect(colour = "black", size = 1),
          axis.title           = element_text(face = "bold"),
          axis.text.x          = element_text(angle = -90, hjust = 0),
          legend.position      = c(1, 1),
          legend.title         = element_blank(),
          legend.justification = c(1, 1)) +
    labs(x = "", y = "READOUTS\n") +
    scale_colour_manual(values = scatter_colour)
}
