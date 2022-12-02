#' Synthetic lethal RNAi screen example data.
#'
#' A dataset containing synthetic lethal RNAi screen data to show how functions work.
#' The variables are as follows (all are character except READOUT):
#'
#' \itemize{
#'   \item PLATE. plate names.
#'   \item MASTER_PLATE. master plate names.
#'   \item WELL_CONTENT_NAME. siRNA targets of wells.
#'   \item EXPERIMENT_TYPE. sample, negative/positive controls.
#'   \item EXPERIMENT_MODIFICATION. experiment conditions, "treatment" or "control".
#'   \item ROW_NAME. row names of plates.
#'   \item COL_NAME. column names of plates.
#'   \item READOUT. screen results.
#' }
#'
#' @docType data
#' @keywords datasets
#' @name example_dt
#' @usage data(example_dt)
#' @format A data.table with 4320 rows and 8 variables
#' @return A data.table containing RANi screen data, the READOUT value has no real biological meaning.
NULL
