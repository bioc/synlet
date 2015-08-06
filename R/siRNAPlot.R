#' Plot siRNA data and quality metrics.
#'
#' Plot the normalized RNAi screen data, row data, control signals and Z' factor.
#' @import magrittr
#' @import ggplot2
#' @import grid
#' @importFrom dplyr filter
#' @importFrom reshape2 melt
#' @param gene gene symbol, case sensitive
#' @param dat synthetic lethal RNAi screen data
#' @param controlsiRNA controlsiRNA could be a vector of several siRNA, including postive/negative control
#' @param FILEPATH path to store the figure
#' @param colour colour used in graphs
#' @param zPrimeMed zPrime factor basing on median
#' @param zPrimeMean, zPrime factor basing on mean
#' @param treatment the treatment condition in EXPERIMENT_MODIFICATION
#' @param control the control condition in EXPERIMENT_MODIFICATION
#' @param normMethod could be a PLATE and negative controls
#' @param width width of the plot
#' @param height height of the plot
#' @return Return the ggplot2 objects in a list, which could be plotted individually.
#' @examples
#' zF_mean <- zFactor(exampleDat, negativeCon = "scrambled control si1",
#'   positiveCon = "PLK1 si1")
#' zF_med <- zFactor(exampleDat, negativeCon = "scrambled control si1",
#'   positiveCon = "PLK1 si1", useMean = FALSE)
#' tem.1 <- siRNAPlot("AAK1", exampleDat,
#'   controlsiRNA = c("lipid only", "scrambled control si1"),
#'   FILEPATH = ".",  zPrimeMed = zF_med, zPrimeMean = zF_mean,
#'   treatment = "treatment", control = "control",
#'   normMethod = c("PLATE", "lipid only", "scrambled control si1"))
#' @export
siRNAPlot <- function(gene, dat, controlsiRNA, FILEPATH = ".",
  colour = rainbow(10), zPrimeMed, zPrimeMean,
  treatment, control, normMethod = c("PLATE"),
  width = 15, height = 14) {
  ## controlsiRNA could be a vector of several siRNA, including postive/negative control;
  ## zPrimeMed, zPrime factor basing on median;
  ## zPrimeMean, zPrime factor basing on mean;
  ## normMethod could be a PLATE and negative controls;
  ## Z factor plots, raw data plots and control plots.
  # browser()
  if(class(dat) != "data.frame") stop("incorrect input!")
  message("+++ Processing", gene, "+++\n")
  tem.1 <- paste("^", gene, " ", sep = "")
  if (length(grep(tem.1, dat$WELL_CONTENT_NAME)) == 0) {
    message("---", gene, "not found!---\n")
  } else {
    tem.match <- paste("^", gene, " ", sep = "")
    tem.11 <- dat[grep(tem.match, dat$WELL_CONTENT_NAME), ][, c("PLATE", "MASTER_PLATE", "WELL_CONTENT_NAME", "EXPERIMENT_MODIFICATION", "READOUT")]
    tem.11$PLATE <- factor(as.character(tem.11$PLATE))

    ## raw siRNA barplots.
    tem.rawBarp <- ggplot(tem.11, aes(x = WELL_CONTENT_NAME, y = READOUT, fill = PLATE)) +
      geom_bar(stat = "identity", position = "dodge") +
      facet_grid(~ EXPERIMENT_MODIFICATION) +
      theme_bw() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      scale_fill_manual(values = colour) +
      theme(axis.text.x = element_text(angle = -90, hjust = 0)) +
      labs(x = "", y = "raw readouts\n", title = "Barplot of raw siRNA screen readouts") +
      theme(legend.title=element_blank())

    ## siRNA control boxplots.
    # browser()
    masterP <- as.character(tem.11[1, "MASTER_PLATE"])
    plates <- unique(as.character(tem.11$PLATE))
    tem.21 <- dplyr::filter(dat, MASTER_PLATE == masterP)
    tem.21 <- tem.21[tem.21$WELL_CONTENT_NAME %in% controlsiRNA, ]
    tem.21$PLATE <- factor(as.character(tem.21$PLATE))
    tem.contBoxp <- ggplot(tem.21, aes(WELL_CONTENT_NAME, READOUT, fill = PLATE)) +
      geom_boxplot() +
      facet_grid(~ EXPERIMENT_MODIFICATION) +
      theme_bw() +
      scale_fill_manual(values = colour) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      labs(x = "", y = "raw readouts\n", title = "Boxplot of raw control readouts") +
      theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
      theme(legend.title=element_blank())

    ## Z' plots
  # browser()
    z.med <- zPrimeMed[plates, ]
    z.me <- zPrimeMean[plates, ]
    z.tot <- cbind(z.med[, 2], z.me[, 2])
    rownames(z.tot) <- rownames(z.med)
    colnames(z.tot) <- c("zPrimeF_median", "zPrimeF_mean")
    tem.31 <- reshape2::melt(z.tot)
    tem.31$Var1 <- factor(tem.31$Var1)
    tem.zp <- ggplot(tem.31, aes(x = Var1, y = value, fill = Var1)) +
      geom_bar(stat = "identity", position = "dodge") +
      facet_grid(~ Var2) +
      theme_bw() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      scale_fill_manual(values = colour) +
      labs(x = "", y = "Z' factor\n", title = "Barplot of Z' factor") +
      theme(axis.text.x = element_text(angle = -90, hjust = 0)) +
      theme(legend.position="none")

    master.norm <- sapply(normMethod, function(X) .ff_masterPlateValue(masterP,
      dat, treatment, control, normMethod = X), simplify = FALSE)
    treat.num <- master.norm[[1]]$treat.num
    cont.num <- master.norm[[1]]$cont.num
    master.norm.c <- do.call(rbind,lapply(names(master.norm),
      function(X) data.frame(siRNA = rownames(master.norm[[X]][[1]]),
        norMethod = X, master.norm[[X]][[1]], check.names = FALSE)))
    tem.match <- paste("^", gene, " ", sep = "")

    # browser()
    siRNA.da <- master.norm.c[grep(tem.match, master.norm.c$siRNA), ]
    siRNA.da <- droplevels(siRNA.da)
  # siRNA.da <- data.frame(siRNA = rownames(siRNA.da), siRNA.da, check.names = FALSE)
  # colnames(siRNA.da)[3:8] <- rep(c("treatment", "control"), each = 3)
    colnames(siRNA.da)[3:8] <- paste(rep(c("treament", "control"),
      times = c(treat.num, cont.num)),
      colnames(siRNA.da)[3:8], sep = "___")
    siRNA.da.m <- reshape2::melt(siRNA.da)
    siRNA.da.m <- data.frame(siRNA.da.m,
      experiments = sapply(as.character(siRNA.da.m$variable),
      function(X) unlist(strsplit(X, "___"))[1]), plates = sapply(as.character(siRNA.da.m$variable), function(X) unlist(strsplit(X, "___"))[2]))

    ## box plot
    tem.boxp <- ggplot(siRNA.da.m, aes(siRNA, value, fill = experiments)) +
      geom_boxplot() +
      theme_bw() +
      theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
      geom_hline(yintercept = 1, linetype = "dashed") +
      facet_grid( ~ norMethod) +
      scale_fill_manual(values = colour) +
      labs(x ='\nsiRNAs', y = 'Normalized Value\n') +
      theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
      labs(title = "Boxplot of normalized value", x = "") +
      theme(legend.title=element_blank())

    ## bar plot
    tem.barp <- ggplot(siRNA.da.m, aes(siRNA, value, fill = experiments)) +
      geom_bar(stat = "identity", position = "dodge") +
      theme_bw() +
      theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
      geom_hline(yintercept = 1, linetype = "dashed") +
      facet_grid(plates ~ norMethod) +
      scale_fill_manual(values = colour) +
      labs(x ='\nsiRNAs', y = 'Normalized Value\n') +
      theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
      labs(title = "Barplot of normalized value", x = "") +
      theme(legend.title=element_blank())

    pdf(file = file.path(FILEPATH, paste(gene, "pdf", sep = ".")), width = width, height = height)
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(2, 6)))
    grid::grid.text("title of this panel", vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1:2))
    print(tem.boxp, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1:3))
    print(tem.barp, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 4:6))
    print(tem.rawBarp, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1:2))
    print(tem.contBoxp, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 3:4))
    print(tem.zp, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 5:6))
    dev.off()

    return(list(tem.boxp = tem.boxp, tem.barp = tem.barp,
      tem.rawBarp = tem.rawBarp, tem.contBoxp = tem.contBoxp,
      tem.zp = tem.zp))

  }
}
