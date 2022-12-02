#' @import data.table
#' @import ggplot2
#' @import patchwork
#' @importFrom magrittr %>% set_rownames set_colnames set_names set_rownames extract
NULL

#- deal with . in magrittr
if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))
