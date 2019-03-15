
# TODO: Add comment
# 
# Author: Vahid Nassiri
###############################################################################

## the dose-response data, we assume first column as dose and secodn as response.
#' simulation-based evalutation of a dose-response model
#' 
#' a function to simulate data based different dose-response model using parameters estimated from a provided pilot study. 
#' The function then simulate data from the estimated model for the given dose levels and number of replications per dose.
#' Some criteria will be compauted which then can be used to compare different settings.
#' 
#' @param pilotData a dataset presenting dose-response data from a pilot study. The first column should give the doses 
#' and the second one should give the response values.
#' @param doseLevels the dose levels which should be used in the simulation study.
#' @param numReplications number of replications for each of thes dose levels for the simulated data.
#' @param numSim number of times that the simulation study should be replicated.
#' @param standardDeviation standard deviation of the generated response.
#' @param EDp scalar in (0,1), indicatign with EDp should be computed to compare different models, default is 0.5 (ED50).
#' @param funcList string vector with models for data generation and fitting, should be selected from 
#' c("linear", "linlog", "exponential", "emax", "sigEmax", "logistic", "betaMod","quadratic").
#' @return a list with the following elements
#' 
#' estEDp a list of length of funcList providing the estimated EDp from models fitted to data generated from each model in funcList
#' realEDp a vector of length funcList, the EDp's computed based on the estimated parameters from different models fitted to pilotData 
#' bestModel a list of length funcList, a frequency table of best selected model for data generated from each model in funcList
#' meanEDp a matrix showing mean of estimated EDp's averaged over numSim replications.
#' biasEDp a matrix showing bias of estimated EDp's averaged over numSim replications.
#' mseEDp a matrix showing MSE of estimated EDp's averaged over numSim replications.
#' varEDp a matrix showing variance of estimated EDp's averaged over numSim replications.
#' relativeBiasEDp a matrix showing relative bias of estimated EDp's averaged over numSim replications.
#' absBiasEDp a matrix showing absolute bias of estimated EDp's averaged over numSim replications.
#' absRelativeBiasEDp a matrix showing absolute bias of estimated EDp's averaged over numSim replications.
#' averagaedAIC a matrix showing AIC's of different models averaged over numSim replications. 
#' quantity2Plot which if needed will be passed to plot method.
#' 
#' The output of simulEvalDRM can be passed to the function plotSimulDRM
#' to plot a heatmap for the desired the quantity of interest. Possible quantities
#' are ("mean", "bias", "mse", "variance", "relativeBias", 
#' 			"absBias", "absRelativeBias")
#' 
#' @importFrom MCPMod genDFdata 
#' @import DoseFinding
#' @examples
#' ## gnerating data, a sample of size 20
#' set.seed(11)
#' doses2Use <-  c(0, 5, 20)
#' numRep2Use <- c(3, 3, 3)
#' generatedData <- cbind(rep(1,sum(numRep2Use)),
#' MCPMod::genDFdata("logistic",c(5, 3, 10, 0.05), doses2Use,
#'  numRep2Use, 1), 
#' 		matrix(rnorm(1*sum(numRep2Use)), sum(numRep2Use), 1))
#' colnames(generatedData) <- c("ID", "dose", "response", "x1")
#' for (iGen in 2:20){
#' genData0 <- cbind(rep(iGen,sum(numRep2Use)), 
#' MCPMod::genDFdata("logistic",c(5, 3, 10, 0.05), doses2Use, 
#' numRep2Use, 1), 
#' 			matrix(rnorm(1*sum(numRep2Use)), sum(numRep2Use), 1))
#' 	colnames(genData0) <- c("ID", "dose", "response", "x1")
#' 	generatedData <- rbind(generatedData, genData0)
#' }
#' simRes <- simulEvalDRM (pilotData = 
#' generatedData[generatedData$ID == 2, c(2,3)], 
#' doseLevels = c(0, 4, 20), 
#' 		numReplications = c(6, 3, 3), numSim = 10, 
#' standardDeviation = 1, EDp = 0.5,
#' 		funcList = c("linlog", "emax", "sigEmax", "logistic"))
#' 
#' @author Vahid Nassiri and Yimer Wasihun.
#' @export
simulEvalDRM <- function(pilotData, doseLevels, numReplications, numSim, standardDeviation, EDp = 0.5,
		funcList = c("linear", "linlog", "exponential", "emax", "sigEmax", "logistic", "betaMod","quadratic")){
	funcList <- match.arg(funcList, c("linear", "linlog", "exponential", "emax", 
					"sigEmax", "logistic", "betaMod","quadratic"),several.ok = TRUE)
	funcListAll <- c("linear", "linlog", "exponential", "emax", "sigEmax", "logistic", "betaMod","quadratic")
	maxDose <- max(doseLevels)
	boundsListAll = list(NULL, NULL, c(0.01*maxDose, maxDose), c(0.01*maxDose, maxDose),
			rbind(c(0.01*maxDose, maxDose), c(0.5, 10)),
			rbind(c(0.01*maxDose, maxDose), c(0.01*maxDose, 0.5*maxDose)), matrix(c(0.05,0.05,4,4), 2), NULL)
	boundsList <- boundsListAll[match(funcList,funcListAll)]
	if (length(doseLevels) != length(numReplications)){
		stop("The length of doseLevels should be equal to length of numReplications")
	}
	if (numSim <0){
		stop("numSim should be a positive number")
	}
	if (floor(numSim) != numSim){
		numSim <- floor(numSim)
		warning(paste("numSim should be an integer, we are using ", numSim))
	}
	## different measures 
	testRes <- NULL
	estEDp <- NULL
	realEDp <- rep(NA, length(funcList))
	meanEDp <- matrix(NA, length(funcList), length(funcList) + 2)
	biasEDp <- matrix(NA, length(funcList), length(funcList) + 2)
	mseEDp <- matrix(NA, length(funcList), length(funcList) + 2)
	varEDp <- matrix(NA, length(funcList), length(funcList) + 2)
	relativeBiasEDp <- matrix(NA, length(funcList), length(funcList) + 2)
	absBiasEDp <- matrix(NA, length(funcList), length(funcList) + 2)
	absRelativeBiasEDp <- matrix(NA, length(funcList), length(funcList) + 2)
	averagaedAIC <- matrix(NA, length(funcList), length(funcList))
	
	for (iTest in 1:length(funcList)){
		## Fit the model
		modelToFit <- funcList[iTest]
		fittedModel <- DoseFinding::fitMod(pilotData[,1], pilotData[,2], data = NULL, 
				model = modelToFit, S = NULL, type = "normal",
				addCovars = ~1, placAdj = FALSE, bnds = boundsList[[which(funcList == modelToFit)]], df = NULL,
				start = NULL, na.action = na.fail, control = NULL,
				addArgs = NULL)
		
		realEDp[iTest] <- localCalcED(model = modelToFit, pars = localSepCoef(fittedModel)$DRpars, 
				p = EDp, maxD = max(doseLevels), EDtype = "continuous", 
				doses = rep(doseLevels, numReplications), off = attr(fittedModel, "off"), 
				scal = attr(fittedModel, "scal"), nodes = attr(fittedModel, "nodes"))
		selectedModel <- rep(NA, numSim)
		## step: generate data from the model
		generatedData <- NULL
		estimatedEDp <- matrix(NA, numSim, length(funcList) + 2)
		modelAIC <- matrix(NA, numSim, length(funcList))
		for (iRep in 1:numSim){
			if (modelToFit == "betaMod"){
				genData <- MCPMod::genDFdata(modelToFit,c(fittedModel$coefs, scal = attr(fittedModel, "scal")), 
						doseLevels, numReplications, standardDeviation)
			}	else{
				genData <- MCPMod::genDFdata(modelToFit,fittedModel$coefs, doseLevels, numReplications, standardDeviation)
			}
			fittedSim <- fitDRMperSubject (subjectData = genData, doseID = 1, responseID = 2, addCovars = ~1,
					funcList = funcList, bounds = boundsList, EDp = EDp, addCovarsVar = FALSE)
			estimatedEDp[iRep, ] <- fittedSim$estEDp
			modelAIC[iRep, ] <- fittedSim$estAIC
			generatedData[[iRep]] <- genData
			selectedModel[iRep] <- funcList[which.min(fittedSim$estAIC)]
		}
		colnames(estimatedEDp) <- c(funcList, "minAIC", "modelAveraging")
		estEDp[[iTest]] <- estimatedEDp 
		testRes [[iTest]] <- table(selectedModel)
		meanEDp[iTest, ] <- apply(estimatedEDp, 2, mean)
		biasEDp[iTest, ] <- apply(estimatedEDp - realEDp[iTest], 2 , mean)
		mseEDp[iTest, ] <- apply((estimatedEDp - realEDp[iTest])^2, 2 , mean)
		varEDp[iTest, ] <- apply(estimatedEDp - realEDp[iTest], 2 , var)
		relativeBiasEDp[iTest, ]  <- apply((estimatedEDp - realEDp[iTest])/realEDp[iTest], 2 , mean)
		absBiasEDp[iTest, ] <- apply(abs(estimatedEDp - realEDp[iTest]), 2 , mean)
		absRelativeBiasEDp[iTest, ] <- apply(abs((estimatedEDp - realEDp[iTest])/realEDp[iTest]), 2 , mean)
		averagaedAIC[iTest, ] <- apply(modelAIC, 2, mean)
	}
	names(estEDp) <- names(testRes) <- funcList
	colnames(meanEDp) <- colnames(biasEDp) <- colnames(mseEDp) <- colnames(varEDp) <- 
			colnames(relativeBiasEDp) <- colnames(absBiasEDp) <- colnames(absRelativeBiasEDp) <- c(funcList, "minAIC", "modelAveraging")
	rownames(meanEDp) <- rownames(biasEDp) <- rownames(mseEDp) <- rownames(varEDp) <- 
			rownames(relativeBiasEDp) <- rownames(absBiasEDp) <- rownames(absRelativeBiasEDp) <- funcList
	rownames(averagaedAIC) <- colnames(averagaedAIC) <- funcList
	toReturn <- list(estEDp = estEDp, realEDp = realEDp, bestModel = testRes, meanEDpAveraged = meanEDp, biasEDpAveraged = biasEDp, 
			mseEDpAveraged = mseEDp, varEDpAveraged = varEDp, relativeBiasEDpAveraged = relativeBiasEDp, 
			absBiasEDpAveraged = absBiasEDp, absRelativeBiasEDpAveraged = absRelativeBiasEDp, avergaedAIC = averagaedAIC)
	class(toReturn) <- "simulDRM"
	return(toReturn)
}
