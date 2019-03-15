# TODO: Add comment
# This function is obtained from sepCoef function in DoseFinding package 
# (https://cran.r-project.org/web/packages/DoseFinding/)
# Author: Vahid Nassiri
###############################################################################


#' a local version of SepCoef in DoseFinding package
#' @param object 
#' @return 
#' 
#' @author Vahid Nassiri (taken from DoseFinding package)
#' @noRd

localSepCoef <- function(object){
	model <- attr(object, "model")
	if(attr(object, "type") == "general")
		return(list(DRpars=object$coefs, covarPars = numeric(0)))
	if(attr(object, "type") == "normal" & object$addCovars == ~1)
		return(list(DRpars=object$coefs, covarPars = numeric(0)))
	## determine the number of parameters (not counting e0 and eMax)
	if(model %in% c("linear","linlog"))
		dim <- 2
	if(model %in% c("quadratic", "exponential", "emax"))
		dim <- 3
	if(model %in% c("sigEmax", "logistic", "betaMod"))
		dim <- 4
	if(model == "linInt")
		dim <- length(attr(object, "nodes"))
	cf <- object$coefs
	p <- length(cf)
	## extract coefficients
	DRpars <- cf[1:dim] # coefs of DR model
	covarPars <- cf[(dim+1):p]
	return(list(DRpars=DRpars, covarPars=covarPars))
}
