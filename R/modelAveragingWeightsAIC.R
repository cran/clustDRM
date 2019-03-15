# TODO: Add comment
# 
# Author: Vahid Nassiri
###############################################################################


#' computes model averaging AIC weights
#' @param estAIC 
#' @return 
#' 
#' @author Vahid Nassiri
#' @noRd

modelAveragingWeightsAIC <- function(estAIC){
	## computes the weights for model averaging
	AICweights0 <- exp(-0.5*(estAIC-min(estAIC)))
	AICweights <- AICweights0/ sum(AICweights0)
	return(AICweights)
}
