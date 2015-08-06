#' Heatmap of all plates
#'
#' Put all individual plates in one graph, values are the readout in experiments.
#' @import magrittr
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom RColorBrewer brewer.pal
#' @param dat synthetic lethal RNAi screen data
#' @param baseSize basic font size used for x/y axis and title for heatmaps
#' @return a ggplot object
#' @examples
#' tem.1 <- plateHeatmap(exampleDat)
#' ggsave("platesHeatmap.pdf", plot = tem.1, width = 500, height = 500, limitsize = FALSE)
#' @export
plateHeatmap <- function(dat, baseSize = 12) {
# browser()
  if(class(dat) != "data.frame") stop("incorrect input!")
  tem.1 <- paste(dat$PLATE, dat$EXPERIMENT_MODIFICATION, sep = " ")
  tem.2 <- cbind(names = tem.1, dat)
  tem.3 <- reshape2::melt(tem.2, id.vars = c("ROW_NAME", "COL_NAME", "names"),
    measure.vars = "READOUT" )
  tem.3$ROW_NAME <- factor(tem.3$ROW_NAME, levels = rev(sort(unique(tem.3$ROW_NAME))))

  myPalette <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, "Spectral")))
  base_size <- baseSize

  ggplot(tem.3, aes(COL_NAME, ROW_NAME)) +
    geom_tile(aes(fill = value), colour = "white") +
    scale_fill_gradientn(colours = myPalette(100)) +
    theme_bw() +
    facet_wrap(~ names) +
    theme(axis.text.x = element_text(size = base_size * 0.5 ), 
      axis.text.y = element_text(size = base_size * 0.5 ), 
      plot.title = element_text(size = base_size * 0.8)) +
    labs(x = "", y = "") +
    coord_equal() +
    scale_x_continuous(expand = c(0, 0)) +
    theme(panel.background=element_rect(fill = "white", colour = "white"))
}
