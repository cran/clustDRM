# TODO: Add comment
# 
# Author: Vahid Nassiri
###############################################################################
#' fits dose-response models to one subject
#' @param subjectData 
#' @param doseID 
#' @param responseID 
#' @param addCovars 
#' @param funcList 
#' @param bounds 
#' @param EDp 
#' @param addCovarsVar 
#' @return
#' 
#' @import DoseFinding 
#' 
#' @author Vahid Nassiri
#' @noRd
fitDRMperSubject <- function(subjectData, doseID, responseID, addCovars,
		funcList, bounds, EDp, addCovarsVar ){
	## For the moment only works for "normal" and "continuous", we should check if the clustering
	## works for the discerete non-normal data, then we can add this option here as well.
	names(subjectData)[doseID] <- "dose"
	names(subjectData)[responseID] <- "response"
#	drcPattern <- match.arg(drcPattern, c("decreasing", "increasing", "maxMin", "minMax"))
#	maxDose <- max(dose)
	## defining the type of functions to fit and boundries based on the pattern of the dose-response curve
	fittedModels <-NULL
	extraCovariates <- NULL
	subjectData <- as.data.frame(subjectData)
	#subjectData[,doseID] <- as.numeric(subjectData[,doseID])
	estEDp <- rep(0, length(funcList))
	estAIC <- rep(0, length(funcList))
	extraCovariatesNames <- all.vars(as.formula(addCovars))
	if (addCovars != as.formula(~1)){
		paramEst <- matrix(NA, length(funcList), length(extraCovariatesNames))
		paramVar <- matrix(NA, length(funcList), length(extraCovariatesNames))
	}

	for (iModel in 1:length(funcList)){
		tmpModel <- DoseFinding::fitMod(dose = "dose", resp = "response", 
				data = subjectData, model = funcList[iModel], S = NULL,
				type = "normal", addCovars = addCovars, bnds = bounds[[iModel]], df = NULL,
				start = NULL)
		fittedModels[[iModel]] <- tmpModel
		estEDp[iModel] <- DoseFinding::ED(tmpModel, p = EDp, EDtype = "continuous")
		estAIC[iModel] <- AIC(tmpModel)
		if (addCovars != as.formula(~1)){
			## store estimates and covariance matrices for extra covariates
			idxAddedCovariate <- which(names(coefficients(tmpModel)) %in% extraCovariatesNames) 
			paramEst[iModel,] <- coefficients(tmpModel)[idxAddedCovariate]
			if (addCovarsVar){
				computedVariance <- try(vcov(tmpModel))
				if (!is.character(computedVariance)){
					paramVar[iModel,] <- diag(computedVariance)[idxAddedCovariate]
				}else{
					addCovarsVar <- FALSE
				}
				
			}
		}
	}
	## estimating EDp based on minimum AIC
	minAICEDp <- estEDp[which.min(estAIC)] 
	## computing model averaging weights based on the AIC
	weightsAIC <- modelAveragingWeightsAIC(estAIC)
	modelAveragingEDp <- sum(weightsAIC * estEDp)
	if (addCovars != as.formula(~1)){
		## combining estimates and standard errors for extra variables
		## formulas obtaiend from Buckland et al. (1997).
		addedCovarsEst <- apply(weightsAIC * paramEst, 2, mean)
		addedCovarsVar <- rep(NA, length(addedCovarsEst))
		if (addCovarsVar){
			paramDiff <- paramEst - matrix(rep(addedCovarsEst,length(funcList)), length(funcList), length(addedCovarsEst),
					byrow = TRUE)
			addedCovarsVar <- (apply(weightsAIC *sqrt(paramVar + (paramDiff^2)), 2, sum))^2
		}
		addedCovariatesModelAveraging <- cbind(addedCovarsEst, sqrt(addedCovarsVar))
		colnames(addedCovariatesModelAveraging) <- c("Est.", "StdErr.")
		addedCovariatesMinAIC <- cbind(paramEst[which.min(estAIC),], sqrt(paramVar[which.min(estAIC),]))
		colnames(addedCovariatesMinAIC) <- c("Est.", "StdErr.")
		extraCovariates <- list(modelAveraging = addedCovariatesModelAveraging, 
				minAIC = addedCovariatesMinAIC)
		class(extraCovariates) <- "drmExtraCovariates"
	}
	## making output
	estEDpFinal <- c(estEDp, minAICEDp, modelAveragingEDp)
	names(estEDpFinal) <- c(funcList, "minAIC", "modelAveragingAIC")
	names(estAIC) <- funcList
	return(list(estEDp = estEDpFinal, fittedModels = fittedModels, estAIC = estAIC, extraCovariates = extraCovariates))
}

