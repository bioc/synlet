#' Plot siRNA data and quality metrics.
#'
#' Plot the normalized RNAi screen data, row data, control signals and Z' factor.
#'
#' @param gene gene symbol, case sensitive
#' @param dta synthetic lethal RNAi screen data
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
#' zF_mean <- zFactor(example_dt, negativeCon = "scrambled control si1", positiveCon = "PLK1 si1")
#' zF_med  <- zFactor(example_dt, negativeCon = "scrambled control si1", positiveCon = "PLK1 si1",
#'                    useMean = FALSE)
#' p01 <- siRNAPlot("AAK1", example_dt,
#'                  controlsiRNA = c("lipid only", "scrambled control si1"),
#'                  FILEPATH = ".",  zPrimeMed = zF_med, zPrimeMean = zF_mean,
#'                  treatment = "treatment", control = "control",
#'                  normMethod = c("PLATE", "lipid only", "scrambled control si1"))
#' @export
#- Need to generate a plot.
siRNAPlot <- function(gene,
    dta,
    controlsiRNA,
    FILEPATH = ".",
    colour   = rainbow(10),
    zPrimeMed,
    zPrimeMean,
    treatment,
    control,
    normMethod = c("PLATE"),
    width      = 15,
    height     = 14) {
  ## controlsiRNA could be a vector of several siRNA, including postive/negative control;
  ## zPrimeMed, zPrime factor basing on median;
  ## zPrimeMean, zPrime factor basing on mean;
  ## normMethod could be a PLATE and negative controls;
  ## Z factor plots, raw data plots and control plots.
  WELL_CONTENT_NAME <- PLATE <- READOUT <- MASTER_PLATE <- rn <- value <- variable <- siRNA <- experiments <- NULL

  message("(==) Processing: ", gene, "\n")

  ex_match <- paste("^", gene, " ", sep = "")

  if (length(grep(ex_match, dta$WELL_CONTENT_NAME)) == 0) {
    message("(II) ", gene, " not found!\n")
  } else {
    dta_2 <- dta[grep(ex_match, WELL_CONTENT_NAME), c("PLATE", "MASTER_PLATE", "WELL_CONTENT_NAME", "EXPERIMENT_MODIFICATION", "READOUT")] %>%
      .[, PLATE := factor(PLATE)]

    #- Raw siRNA barplots.
    p_raw_barpl <- ggplot(as.data.frame(dta_2), aes( WELL_CONTENT_NAME, READOUT, fill = PLATE)) +
      geom_bar(stat = "identity", position = "dodge") +
      facet_grid(~ EXPERIMENT_MODIFICATION) +
      theme_bw(14) +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_rect(colour = "black", size = 1),
            axis.text.x      = element_text(angle = -90, hjust = 0),
            legend.title     = element_blank()) +
      scale_fill_manual(values = colour) +
      labs(x = "", y = "raw readouts", title = "Barplot of raw screen readouts")

    dta_3 <- dta[MASTER_PLATE == unique(dta_2$MASTER_PLATE)] %>%
      .[WELL_CONTENT_NAME %in% controlsiRNA] %>%
      .[, PLATE := factor(PLATE)]

    p_cont_boxplot <- ggplot(as.data.frame(dta_3), aes(WELL_CONTENT_NAME, READOUT, fill = PLATE)) +
      geom_boxplot() +
      facet_grid(~ EXPERIMENT_MODIFICATION) +
      theme_bw(14) +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_rect(colour = "black", size = 1),
            axis.text.x      = element_text(angle = -45, hjust = 0),
            legend.title     = element_blank()) +
      scale_fill_manual(values = colour) +
      labs(x = "", y = "raw readouts", title = "Boxplot of raw control readouts")

    z_med  <- zPrimeMed[as.character(unique(dta_2$PLATE)), ]
    z_mean <- zPrimeMean[as.character(unique(dta_2$PLATE)), ]
    z_comb <- cbind(z_med[, 2], z_mean[, 2]) %>%
      set_rownames(rownames(z_med)) %>%
      set_colnames(c("zPrimeF_median", "zPrimeF_mean")) %>%
      as.data.table(keep.rownames = TRUE) %>%
      #- Set the var.
      melt(measure.vars = c("zPrimeF_median", "zPrimeF_mean"))

    #- Z primer plots.
    p_zprimer_barpl <- ggplot(as.data.frame(z_comb), aes(rn, value, fill = rn)) +
      geom_bar(stat = "identity", position = "dodge") +
      facet_grid(~ variable) +
      theme_bw(14) +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_rect(colour = "black", size = 1),
            axis.text.x      = element_text(angle = -90, hjust = 0),
            legend.title     = element_blank()) +
      scale_fill_manual(values = colour) +
      labs(x = "", y = "Z' factor", title = "Barplot of Z' factor")

    master_norm <- sapply(normMethod,
                          function(X) .ff_masterPlateValue(unique(dta_2$MASTER_PLATE),
                                                           dta,
                                                           treatment,
                                                           control,
                                                           normMethod = X), simplify = FALSE)
    n_treat_plate <- master_norm[[1]]$n_treat_plate
    n_cont_plate  <- master_norm[[1]]$n_cont_plate
    master_norm_c <- do.call(rbind,lapply(names(master_norm),
                                          function(x) data.frame(siRNA       = rownames(master_norm[[x]][[1]]),
                                                                 norMethod   = x,
                                                                 master_norm[[x]][[1]],
                                                                 check.names = FALSE)))

    dta_4 <- master_norm_c[grep(ex_match, master_norm_c$siRNA), ] %>%
      droplevels %>%
      as.data.table %>%
      setnames(2 + seq_len(n_treat_plate + n_cont_plate), paste(rep(c("treament", "control"), times = c(n_treat_plate, n_cont_plate)),
                                                                names(.)[2 + seq_len(n_treat_plate + n_cont_plate)],
                                                                sep = "__")) %>%
      melt(measure.vars = names(.)[-c(1:2)]) %>%
      .[, c("experiments", "plates") := tstrsplit(variable, split = "__")]

    #- box plot
    p_norm_boxpl <- ggplot(as.data.frame(dta_4), aes(siRNA, value, fill = experiments)) +
      geom_boxplot() +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text.x      = element_text(angle = -45, hjust = 0),
            legend.title     = element_blank(),
            panel.background = element_rect(colour = "black", size = 1)) +
      geom_hline(yintercept = 1, linetype = "dashed") +
      facet_grid( ~ norMethod) +
      scale_fill_manual(values = colour) +
      labs(x ="", y = "Normalized Value", title = "Boxplot of normalized value")

    #- bar plot
    p_norm_barpl <- ggplot(as.data.frame(dta_4), aes(siRNA, value, fill = experiments)) +
      geom_bar(stat = "identity", position = "dodge") +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text.x      = element_text(angle = -45, hjust = 0),
            legend.title     = element_blank(),
            panel.background = element_rect(colour = "black", size = 1)) +
      geom_hline(yintercept = 1, linetype = "dashed") +
      facet_grid(plates ~ norMethod) +
      scale_fill_manual(values = colour) +
      labs(x ="", y = "Normalized Value", title = "Barplot of normalized value")

    q01 <- (p_norm_boxpl | p_norm_barpl) / (p_raw_barpl | p_cont_boxplot | p_zprimer_barpl)

    ggsave(file.path(FILEPATH, paste(gene, "pdf", sep = ".")), q01, width = width, height = height)

    return(list(p_norm_boxpl    = p_norm_boxpl,
                p_norm_barpl    = p_norm_barpl,
                p_raw_barpl     = p_raw_barpl,
                p_cont_boxplot  = p_cont_boxplot,
                p_zprimer_barpl = p_zprimer_barpl))

  }
}
