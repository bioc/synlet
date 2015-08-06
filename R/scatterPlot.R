#' Scatter plot of RNAi screen results
#'
#' Produce a single plot for readous of each plate, with the option of highlighting
#' specific signals, like positive/negative controls.
#'
#' @import magrittr
#' @import ggplot2
#' @param dat synthetic lethal RNAi screen data
#' @param controlOnly whether or not to plot control wells only
#' @param colour colour for different signals
#' @param ... positive/negative signals, must be specified
#' @return a ggplot object
#' @examples
#' scatterPlot(exampleDat, controlOnly = FALSE, colour = rainbow(10),
#'   "PLK1 si1", "scrambled control si1", "lipid only")
#' @export
scatterPlot <- function(dat, controlOnly = FALSE, colour, ...) {
  if(class(dat) != "data.frame") stop("incorrect input!")
  xx.1 <- as.character(dat[, "WELL_CONTENT_NAME"])
  ## specific signals
  if (controlOnly == FALSE) {
    xx.1[grep(paste(c(...), collapse="|"), xx.1, invert=TRUE)] <- "knock-down signals"
    ## location of controll wells;
    data.f2 <- cbind(dat, condition = xx.1)
    } else {
    xx.control <- grep(paste(c(...), collapse="|"), xx.1)
    xx.2 <- xx.1[xx.control]
    data.f2 <- cbind(dat[xx.control, ], condition = xx.2)
}

  data.f2 <- data.f2[order(data.f2$EXPERIMENT_MODIFICATION, decreasing = TRUE),]
  data.f2$PLATE <- factor(data.f2$PLATE)
  xx.1 <- paste(data.f2$PLATE, data.f2$EXPERIMENT_MODIFICATION, sep = " ") %>%
    factor
  data.f2 <- cbind(data.f2, names = xx.1)
  ggplot(data.f2, aes(names, READOUT, colour = condition)) +
    geom_point() +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(axis.title = element_text(face = "bold")) +
    theme(axis.text.x = element_text(angle = -90, hjust = 0)) +
    labs(y = "READOUTS\n") +
    theme(legend.position = c(1,1), legend.justification = c(1,1)) +
    scale_colour_manual(values = colour)
}
