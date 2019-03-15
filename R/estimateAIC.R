# TODO: Add comment
# 
# Author: Vahid Nassiri
###############################################################################
## Note that if we want to compute AICc, then:
## logLikelihoodDFc <- logLikelihoodDF + (2*logLikelihoodDF*(logLikelihoodDF+1))/(n-logLikelihoodDF-1)
## and then this logLikelihoodDFc should be used for correction.

#' estimates AIC for a given model, estimated parameters, and set of data
#' @param generatedParams 
#' @param fittedModel 
#' @return 
#' 
#' @author Vahid Nassiri
#' @noRd

estimateAIC <- function(generatedParams, fittedModel){
	#parEst <- coef(fittedModel)
	inputData <- attr(fittedModel, "data")
	n <- nrow(inputData)
	modelToUse <- attr(fittedModel, "model")
	scal <- attr(fittedModel, "scal")
	off <- attr(fittedModel, "off")
	nodes <- attr(fittedModel, "nodes")
	## Obtain the part related to the added covars from the data 
	addedCovars <- attr(fittedModel, "addCovars")
	yhatAddedCovars <- 0
	if (addedCovars != as.formula(~1)){
		numAddCovars <- length(all.vars(as.formula(addedCovars)))
		addevCovarsData <- matrix(unlist(inputData[,3:(2+numAddCovars)]), n, numAddCovars)
		yhatAddedCovars <- addevCovarsData%*%generatedParams[names(generatedParams) %in% all.vars(as.formula(addedCovars))] 
	}
	if (modelToUse == "linear"){
		yhat <- generatedParams[1] + generatedParams[2] * inputData[,1] + yhatAddedCovars
	}
	if (modelToUse == "linlog"){
		yhat <- generatedParams[1] + (generatedParams[2] *log(inputData[,1] + off))+ yhatAddedCovars
	}
	if (modelToUse == "exponential"){
		yhat <- generatedParams[1]+ (generatedParams[2] * (exp(inputData[,1]/ generatedParams[3])-1)) + yhatAddedCovars
	}
	if (modelToUse == "emax"){
		yhat <- generatedParams[1] + (generatedParams[2] * (inputData[,1]/(generatedParams[3] + inputData[,1]))) + yhatAddedCovars
	}
	if (modelToUse == "sigEmax"){
		yhat <- generatedParams[1] + (generatedParams[2] * ((inputData[,1]^generatedParams[4])/((generatedParams[3]^generatedParams[4]) + 
							(inputData[,1]^generatedParams[4])))) + yhatAddedCovars
	}
	if (modelToUse == "logistic"){
		yhat <- generatedParams[1] + (generatedParams[2]/(1+ exp((generatedParams[3] - inputData[,1])/generatedParams[4]))) + yhatAddedCovars
	}
	if (modelToUse == "quadratic"){
		yhat <- generatedParams[1] + (generatedParams[2] * inputData[,1]) + (generatedParams[3] * (inputData[,1]^2)) + yhatAddedCovars
	}
	if (modelToUse == "betaMod"){
		B <- ((generatedParams[3] + generatedParams[4])^(generatedParams[3] + generatedParams[4])) / 
				((generatedParams[3]^generatedParams[3])*(generatedParams[4]^generatedParams[4]))
		yhat <- generatedParams[1] + 
				(generatedParams[2] * B * ((inputData[,1]/scal)^generatedParams[3]) * ((1-(inputData[,1]/scal))^generatedParams[4])) + 
				yhatAddedCovars
	}
	## computing residual sum of squares of the normal regression model
	RSS <- sum((inputData[,2] - yhat)^2)
	## estimating error variance of the normal regression model
	## Note that, sig2 is based on yhat, and yhat is based on generated parameters,
	## so even estimated error variance is different every time.
	sig2 <- RSS/n
	logLikelihood <- -n/2*(log(2*pi) + 1 + log(sig2))
	## the +1 is because we have estimated an extra parameter: sigma2
	logLikelihoodDF <- length(generatedParams)+1 
	estAIC <-   -2*as.vector(logLikelihood) + 2*(logLikelihoodDF)
	toReturn <- c(estAIC, logLikelihood, logLikelihoodDF)
	names(toReturn) <- c("AIC", "logLik", "logLikDF")
	return(toReturn)
}
	 