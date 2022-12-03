#' @import data.table
#' @import ggplot2
#' @import patchwork
#' @importFrom magrittr %>% set_rownames set_colnames set_names set_rownames extract inset
#' @importFrom grDevices colorRampPalette rainbow
#' @importFrom methods is
#' @importFrom utils write.table
#' @importFrom stats mad median medpolish na.omit p.adjust phyper sd t.test
NULL

#- deal with . in magrittr
if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))
