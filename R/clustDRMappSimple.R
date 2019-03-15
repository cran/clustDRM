# TODO: Add comment
# 
# Author: Vahid Nassiri
###############################################################################



#' launch the simpler version of the shiny app for an easier use of the package 
#' 
#' @details a simpler version of the shiny app made for an easy use of the functionalities 
#' of the clustDRM package. It can be launched using command clustDRMappSimple(). 
#' It imports the data (in csv format),
#' performs the clustering on it (for monotone patterns using E2
#' test, and for general patterns using ORICC2 and MCT). It also
#' can estimate EDp for different dose-response curves using 
#' appropriate models. Plotting dose-response curves are also
#' possible for any of these operations. 
#' 
#' @import shiny 
#' @importFrom parallel detectCores
#' @author Vahid Nassiri and Yimer Wasihun.
#' @export
clustDRMappSimple <- function() {
	options(shiny.maxRequestSize = 30*1024^2)
	shiny::runApp( appDir = system.file("uiSimpler", package = "clustDRM"))
}
