# TODO: Add comment
# 
# Author: Vahid Nassiri
###############################################################################
#' fitting dose-response model according to the identified pattern.
#' 
#' function to fit several dose-response candidate models according to the identified pattern, and combine their 
#' results using model selection and/or model averaging.
#' 
#' @param inputDataset a data frame containing the input dataset, it should at least include dose, response, and ID
#' @param dose either a single string or a scalar, indicating the name of the dose column or its index.
#' @param response either a single string or a scalar, indicating the name of the response column or its index.
#' @param ID either a single string or a scalar, indicating the name of the ID column or its index.
#' @param subsettingID a vector of ID's of the subjects, in case one wants to fit the models only to a subset of the data. Default is NULL, 
#' i.e., all the subjects in the inputDataset will be used.
#' @param transform single string indicating what kind of transform should be applied on the response data. 
#' It takes "none" (no transform, dafault), "log" (natural log), "sRoot (square root), and "qRoot" (cubic root), and 
#' "boxcox" (Box-Cox transformation).
#' @param addCovars formula specifying extra linear covariate, e.g., ~x1+x2
#' @param patternClusters a vector of the same length as the number of rows in inputData (number of subjects) indicating a 
#' pattern for each subject. Note that the keywords which are recognized are: "increasing", "decreasing", "flat", "complete", and 
#' "up down max at x" and "down up min at x", which x is one of the doses. The "flat" and "complete" patterns would not be considered.
#' @param EDp scalar in (0,1), indicatign with EDp should be computed, default is 0.5 (ED50).
#' @param addCovarsVar logical variable (TRUE as default), idicating whether the variance of the extra covariates 
#' preesented in  addCovars (unless it is only intercept) should also be computed or not.
#' @param alpha scalar in (0,1), level of significance with default alpha = 0.05.
#' @param na.rm logical variable indicatign whether missing values should be removed (TRUE) or not (FALSE, default)
#' @param imputationMethod signle string taking calues from "mean" (default), and "median", which indicates how the missing values should be
#' treated. "mean" would replace them with the mean of the observed ones, and "median" will use median of them for imputation.
#' @param nCores scalar, indicating the number of cores should be used to perform LRT and MCT tests. Default is 1 which means sequantial 
#' computation (no prallel computation).
#' 
#' @details Note that the dose column of the inputDataset should be a numeric variable.
#' 
#' 
#' @return an object fo class fittedDRM which is a list with the following objects:
#'  fittedModels the outcome of DoseFinding::fitMod for all the suitable models
#'  estAICNonmonotone: the computed AIC for the models fitted to the subjects with a non-monotone pattern
#'  estEDpNonmonotone: the computed EDp for the models fitted to the subjects with a non-monotone pattern
#'  estAICMonmonotone: the computed AIC for the models fitted to the subjects with a monotone pattern
#'  estEDpMonmonotone: the computed EDp for the models fitted to the subjects with a monotone pattern
#'  extraCovarsMonotone: if any extra covariates are added to the model their estimates and possibly standard errors 
#'  (if addCovarsVar = TRUE) are gievn for subjects with monotone pattern.
#' extraCovarsNonmonotone: if any extra covariates are added to the model their estimates and possibly standard errors 
#'  (if addCovarsVar = TRUE) are gievn for subjects with non-monotone pattern.
#' 
#' @examples
#' ## gnerating data
#' set.seed(11)
#' doses2Use <-  c(0, 5, 20)
#' numRep2Use <- c(6, 3, 3)
#' generatedData <- cbind(rep(1,sum(numRep2Use)),
#' MCPMod::genDFdata("logistic",c(5, 3, 10, 0.05), doses2Use, 
#' numRep2Use, 1), 
#' 		matrix(rnorm(1*sum(numRep2Use)), sum(numRep2Use), 1))
#' colnames(generatedData) <- c("ID", "dose", "response", "x1")
#' for (iGen in 2:15){
#' 	genData0 <- cbind(rep(iGen,sum(numRep2Use)),
#' MCPMod::genDFdata("logistic",c(5, 3, 10, 0.05), doses2Use, 
#' numRep2Use, 1), matrix(rnorm(1*sum(numRep2Use)), 
#' sum(numRep2Use), 1))
#' 	colnames(genData0) <- c("ID", "dose", "response", "x1")
#' 	generatedData <- rbind(generatedData, genData0)
#' }
#' ## transforming it for clustering
#' toInput <- inputDataMaker(2, 3, 1, generatedData)
#' ## general pattern clustering
#' generalPatternClust <- generalPatternClustering(
#' inputData = toInput$inputData, colsData = toInput$colsData ,
#' colID = toInput$colID, doseLevels = toInput$doseLevels, 
#' numReplications = toInput$numReplicates, na.rm = FALSE, 
#' imputationMethod = "mean", ORICC = "two", transform = "none",
#' plotFormat = "eps", LRT = TRUE, MCT = TRUE,
#' 		adjustMethod = "BH", nPermute = 100, useSeed = NULL, 
#' theLeastNumberOfMethods = 2, alpha = 0.05, nCores = 1)
#' ## fitDRM 
#' fittedModel <- fitDRM (inputDataset = generatedData, dose = 2, 
#' response = 3, ID = 1, subsettingID = NULL, 
#' 		transform = c("none"), addCovars = ~x1, 
#' 		patternClusters = 
#' generalPatternClust$clusteringORICC2Results$clusteringResultsORICC2, 
#' 		EDp = 0.5, addCovarsVar = TRUE, alpha = 0.05, na.rm = FALSE, 
#' imputationMethod = c("mean"), nCores = 1)
#' 
#' @seealso \href{https://www.rdocumentation.org/packages/DoseFinding/versions/0.9-16/topics/fitMod}{DoseFinding}
#' @import DoseFinding
#' @import doParallel
#' @import foreach
#' @import parallel
#' @importFrom parallel detectCores
#' @importFrom stats AIC aov as.formula coefficients na.fail p.adjust rnorm sd uniroot var vcov
#' @author Vahid Nassiri, and Yimer Wasihun
#' @export

fitDRM <- function(inputDataset, dose, response, ID, subsettingID = NULL, 
		transform = c("none", "log","sRoot", "qRoot", "boxcox"), addCovars = ~1, patternClusters, 
		EDp = 0.5, addCovarsVar = TRUE, alpha = 0.05,
		na.rm = FALSE, imputationMethod = c("mean", "median"), nCores = 1){
	
	i <- iLRT <- iMCT <- NULL
	
## making a list out of the inputDataset with data for each subject as its elements
	
	colID <- ID
	doseColumn <- dose 
	responseCol <- response
	if (is.character(ID)){
		colID <- which(names(inputDataset) == ID)
	}
	
	if (is.character(dose)){
		doseColumn <- which(names(inputDataset) == dose)
	}
	
	if (is.character(response)){
		responseCol <- which(names(inputDataset) == response)
	}
	## treating missing value
	idxMiss <- unique(inputDataset[,colID][which(is.na(inputDataset[,responseCol])== TRUE)])
	if (na.rm){
		inputDataset <- inputDataset[-which((inputDataset[,colID] %in% idxMiss) == TRUE),]
	}else{
		for (iMiss in seq_along(idxMiss)){
			missSubject <- inputDataset[inputDataset[,colID] == idxMiss[iMiss], responseCol]
			inputDataset[inputDataset[,ID] == idxMiss[iMiss] & is.na(inputDataset[,responseCol]), responseCol] <- 
					get(imputationMethod)(missSubject, na.rm = TRUE)
		}
	}
	## transforming data if needed
	if (transform == "none"){
		inputDataset[,responseCol] <- inputDataset[,responseCol]
	}else if(transform == "log"){
		if (any(unlist(inputDataset[,responseCol]) == 0 )){
			inputDataset[,responseCol] <- log(inputDataset[,responseCol] + 1)
			warning("Due to zero responses, log(x+1) is used instead of log.")
		}
		if(any(unlist(inputDataset[,responseCol]) <0)){
			stop("log transofrm is not possible with negative responses.")
		}
		if(all(unlist(inputDataset[,responseCol])>0)){
			inputDataset[,responseCol] <- log(inputDataset[,responseCol])
		}
	}else if (transform == "sRoot"){
		if (any(unlist(inputDataset[,responseCol])<0)){
			stop("square root is not possible with negative responses.")
		}else{
			inputDataset[,responseCol] <- sqrt(inputDataset[,responseCol])
		}
	}else if (transform == "qRoot"){
		inputDataset[,responseCol] <- sign(inputDataset[,responseCol]) * (abs(inputDataset[,responseCol])^(1/3))
	}else if(transform  == "boxcox"){
		inputDataset[,responseCol] <- apply(inputDataset[,responseCol], 2, boxCoxTransform)
	}
	## making a list out of inputDataset 	
	inputDataList <- split(inputDataset, f = unlist(inputDataset[,colID]))
	if (!is.null(subsettingID)){
		## check to see if all the subsettingID list are actually in the data after possibly removing NA's
		if (!all(subsettingID %in% names(inputDataList))){
			stop("Some of ID's in subsettingID are not in the dataset, either they were not there from the beginning or they are removed because na.rm = TRUE")
		}
		inputDataList <- inputDataList[which(names(inputDataList) %in% subsettingID == TRUE)]
	}
	if (length(patternClusters) != length(inputDataList)){
		stop("The length of patternClusters should be equal to number of rows in inputData.")
	}
	#originalInputDataset <- inputDataset
	## finding maximum dose
	maxDose <- max(inputDataset[,doseColumn])
	## providing list of candidatre models if we pattern is monotone or not
	funcListMonotone <- c("linear", "linlog", "exponential", "emax", "sigEmax", "logistic")
	funcListNonmonotone <- c("betaMod","quadratic")
	## using DoseFinding package, providing the boundy values for monotone and non-monotone patterns
	boundsMonotone <- list(NULL, NULL, c(0.01*maxDose, maxDose), c(0.01*maxDose, maxDose),
			rbind(c(0.01*maxDose, maxDose), c(0.5, 10)),
			rbind(c(0.01*maxDose, maxDose), c(0.01*maxDose, 0.5*maxDose)))
	boundsNonmonotone <- list(matrix(c(0.05,0.05,4,4), 2), NULL)
	## finding the idx of monotone and non-monotone (except complete.profile)
	idxMonotone <- which(patternClusters == "decreasing" | patternClusters == "increasing")
	idxNonmonotone <- which(patternClusters != "decreasing" & patternClusters != "increasing"
					& patternClusters != "complete" & patternClusters != "flat")
	## fitting t he model
	#resultsNonmonotoneBoot <- NULL
	#resultsMonotoneBoot <- NULL
	resultsNonmonotone <- NULL
	resultsMonotone <- NULL
	#EDpBoot <- NULL
	############ changhe colID to colsData, just to be sure
	if (nCores > parallel::detectCores()){
		stop(paste("Your system only has ", detectCores, "you cannot specify an nCores larger than that."))
	}else if(nCores == 1){
		if (length(idxMonotone) > 0){
			resultsMonotone <- mapply(fitDRMperSubject, subjectData = inputDataList[idxMonotone], 
					doseID = rep(list(doseColumn), length(idxMonotone)), 
					responseID = rep(list(responseCol), length(idxMonotone)),
					addCovars = rep(list(addCovars), length(idxMonotone)),
					funcList = rep(list(funcListMonotone), length(idxMonotone)), 
					bounds = rep(list(boundsMonotone), length(idxMonotone)), 
					EDp = rep(list(EDp), length(idxMonotone)),
					addCovarsVar = rep(list(addCovarsVar), length(idxMonotone)),
					SIMPLIFY = FALSE)
		}
		if (length(idxNonmonotone) > 0){
			resultsNonmonotone <- mapply(fitDRMperSubject, subjectData = inputDataList[idxNonmonotone], 
					doseID = rep(list(doseColumn), length(idxNonmonotone)), 
					responseID = rep(list(responseCol), length(idxNonmonotone)),
					addCovars = rep(list(addCovars), length(idxNonmonotone)),
					funcList = rep(list(funcListNonmonotone), length(idxNonmonotone)), 
					bounds = rep(list(boundsNonmonotone), length(idxNonmonotone)), 
					EDp = rep(list(EDp), length(idxNonmonotone)),
					addCovarsVar = rep(list(addCovarsVar), length(idxNonmonotone)),
					SIMPLIFY = FALSE)
		}
#		if (bootAveraging){
#			resultsMonotoneBoot <- mapply(EDpBootstrapModelAveraging, 
#					inputDataList[idxMonotone], 
#					doseID = rep(list(doseColumn), length(idxMonotone)), 
#					responseID = rep(list(responseCol), length(idxMonotone)),
#					addCovars = rep(list(addCovars), length(idxMonotone)),
#					funcList = rep(list(funcListMonotone), length(idxMonotone)), 
#					bounds = rep(list(boundsMonotone), length(idxMonotone)), 
#					EDp = rep(list(EDp), length(idxMonotone)), 
#					nBoot = rep(list(nBoot), length(idxMonotone)), 
#					alpha = rep(list(alpha), length(idxMonotone)), 
#					useSeed = rep(list(useSeed), length(idxMonotone)),
#					SIMPLIFY = FALSE)
#			resultsNonmonotoneBoot <- mapply(EDpBootstrapModelAveraging, 
#					inputDataList[idxNonmonotone], 
#					doseID = rep(list(doseColumn), length(idxNonmonotone)), 
#					responseID = rep(list(responseCol), length(idxNonmonotone)),
#					addCovars = rep(list(addCovars), length(idxNonmonotone)),
#					funcList = rep(list(funcListNonmonotone), length(idxNonmonotone)), 
#					bounds = rep(list(boundsNonmonotone), length(idxNonmonotone)), 
#					EDp = rep(list(EDp), length(idxNonmonotone)), 
#					nBoot = rep(list(nBoot), length(idxNonmonotone)), 
#					alpha = rep(list(alpha), length(idxNonmonotone)), 
#					useSeed = rep(list(useSeed), length(idxNonmonotone)),
#					SIMPLIFY = FALSE)
#			names(resultsNonmonotoneBoot) <- names(inputDataList)[idxNonmonotone]
#			#names(resultsMonotone) <- inputData[,colID][idxMonotone]
#			names(resultsMonotoneBoot) <- names(inputDataList)[idxMonotone]
#			class(resultsNonmonotoneBoot) <- "bootstrapModelAveraging"
#			class(resultsMonotoneBoot) <- "bootstrapModelAveraging"
#			## filling in a matrix for estimated EDp's with the CI
#			EDpBoot <- matrix(NA, length(inputDataList), 3)
#			for (iEDp in seq_along(idxMonotone)){
#				EDpBoot[idxMonotone[iEDp],] <- resultsMonotoneBoot[[iEDp]]$bootstrapEDp
#			}
#			for (iEDp in seq_along(idxNonmonotone)){
#				EDpBoot[idxNonmonotone[iEDp],] <- resultsNonmonotoneBoot[[iEDp]]$bootstrapEDp
#			}
#			EDpBoot <- data.frame(names(inputDataList), EDpBoot)
#			colnames(EDpBoot) <- c(colnames(inputDataset)[colID], "median", paste(alpha/2*100,"%", sep = ""), 
#					paste((1-alpha/2) * 100,"%", sep = ""))
#		}
	}else{
		if (Sys.info()['sysname'] == "Windows"){
			cl <- makeCluster(nCores, type="PSOCK", outfile="")  
		}else{
			cl <- makeCluster(nCores, type="FORK", outfile="")  
		}
		registerDoParallel(cl)  
		if (length(idxMonotone) > 0){
			resultsMonotone <- foreach(iMono=seq_along(idxMonotone)) %dopar% 
					fitDRMperSubject(subjectData = inputDataList[[idxMonotone[iMono]]], doseID = doseColumn, responseID = responseCol, 
							addCovars = addCovars, 
							funcList = funcListMonotone, bounds = boundsMonotone, EDp = EDp, addCovarsVar = addCovarsVar)
		}
		if (length(idxNonmonotone) > 0){
			resultsNonmonotone <- foreach(iNonmono=seq_along(idxNonmonotone)) %dopar% 
					fitDRMperSubject(subjectData = inputDataList[[idxNonmonotone[iNonmono]]], doseID = doseColumn, responseID = responseCol, 
							addCovars = addCovars, 
							funcList = funcListNonmonotone, bounds = boundsNonmonotone, EDp = EDp, addCovarsVar = addCovarsVar)
		}
#		if (bootAveraging){
#			resultsMonotoneBoot <- foreach(iMono=seq_along(idxMonotone)) %dopar% 
#					EDpBootstrapModelAveraging(subjectData = inputDataList[[idxMonotone[iMono]]], doseID = doseColumn, responseID = responseCol, 
#							addCovars = addCovars,
#							funcList = funcListMonotone, bounds = boundsMonotone, EDp = EDp, nBoot = nBoot, alpha = alpha, useSeed = useSeed)
#			resultsNonmonotoneBoot <- foreach(iNonmono=seq_along(idxNonmonotone)) %dopar% 
#					EDpBootstrapModelAveraging(subjectData = inputDataList[[idxNonmonotone[iNonmono]]], doseID = doseColumn, responseID = responseCol, 
#							addCovars = addCovars, 
#							funcList = funcListNonmonotone, bounds = boundsNonmonotone, EDp = EDp, nBoot = nBoot, alpha = alpha, useSeed = useSeed)
#		}
		stopCluster(cl)
#		if (bootAveraging){
#			names(resultsNonmonotoneBoot) <- names(inputDataList)[idxNonmonotone]
#			#names(resultsMonotone) <- inputData[,colID][idxMonotone]
#			names(resultsMonotoneBoot) <- names(inputDataList)[idxMonotone]
#			class(resultsNonmonotoneBoot) <- "bootstrapModelAveraging"
#			class(resultsMonotoneBoot) <- "bootstrapModelAveraging"
#			## filling in a matrix for estimated EDp's with the CI
#			EDpBoot <- matrix(NA, length(inputDataList), 3)
#			for (iEDp in seq_along(idxMonotone)){
#				EDpBoot[idxMonotone[iEDp],] <- resultsMonotoneBoot[[iEDp]]$bootstrapEDp
#			}
#			for (iEDp in seq_along(idxNonmonotone)){
#				EDpBoot[idxNonmonotone[iEDp],] <- resultsNonmonotoneBoot[[iEDp]]$bootstrapEDp
#			}
#			EDpBoot <- data.frame(names(inputDataList), EDpBoot)
#			colnames(EDpBoot) <- c(colnames(inputDataset)[colID], "median", paste(alpha/2*100,"%", sep = ""), 
#					paste((1-alpha/2) * 100,"%", sep = ""))
#		}
	}
	## extract results from the outcome
	fittedModelsMonotone <- NULL
	estAICMonotone <- matrix(0, length(idxMonotone), 6)
	estEDpMonotone <- matrix(0, length(idxMonotone), 8)
	fittedModelsNonmonotone <- NULL
	estAICNonmonotone <- matrix(0, length(idxNonmonotone), 2)
	estEDpNonmonotone <- matrix(0, length(idxNonmonotone), 4)
	extraVarsMonotone <- NULL
	extraVarsNonmonotone <- NULL
	
	if (addCovars != as.formula(~1)){
		extraVarsMonotone <- matrix(NA, length(idxMonotone), length(all.vars(as.formula(addCovars))) * 4)
		for (iMono in seq_along(idxMonotone)){
			fittedModelsMonotone[[iMono]] <- resultsMonotone[[iMono]]$fittedModels
			estAICMonotone[iMono,] <- resultsMonotone[[iMono]]$estAIC
			estEDpMonotone[iMono,] <- resultsMonotone[[iMono]]$estEDp
			extraVarsMonotone[iMono, ] <- c(c(resultsMonotone[[iMono]]$extraCovariates$minAIC)
					,c(resultsMonotone[[iMono]]$extraCovariates$modelAveraging))
		}
		colnames(estAICMonotone) <- funcListMonotone
		colnames(estEDpMonotone) <- c(funcListMonotone, "minAIC", "modelAveragingAIC")
		
		## Non-monotone
		extraVarsNonmonotone <- matrix(NA, length(idxNonmonotone), length(all.vars(as.formula(addCovars))) * 4)
		for (iNonmono in seq_along(idxNonmonotone)){
			fittedModelsNonmonotone[[iNonmono]] <- resultsNonmonotone[[iNonmono]]$fittedModels
			estAICNonmonotone[iNonmono,] <- resultsNonmonotone[[iNonmono]]$estAIC
			estEDpNonmonotone[iNonmono,] <- resultsNonmonotone[[iNonmono]]$estEDp
			extraVarsNonmonotone[iNonmono, ] <- c(c(resultsNonmonotone[[iNonmono]]$extraCovariates$minAIC)
					,c(resultsNonmonotone[[iNonmono]]$extraCovariates$modelAveraging))
		}
		colnames(estAICNonmonotone) <- funcListNonmonotone
		colnames(estEDpNonmonotone) <- c(funcListNonmonotone, "minAIC", "modelAveragingAIC")
		## adding ID's
		fittedModelsFinal<- as.list(rep(NA, length(inputDataList)))
		fittedModelsFinal[idxMonotone] <- fittedModelsMonotone
		fittedModelsFinal[idxNonmonotone] <- fittedModelsNonmonotone
		estAICNonmonotone <- data.frame(names(inputDataList)[idxNonmonotone],estAICNonmonotone)
		estEDpNonmonotone <- data.frame(names(inputDataList)[idxNonmonotone],estEDpNonmonotone)
		estAICMonotone <- data.frame(names(inputDataList)[idxMonotone],estAICMonotone)
		estEDpMonotone <- data.frame(names(inputDataList)[idxMonotone],estEDpMonotone)
		extraVarsMonotone <- data.frame(names(inputDataList)[idxMonotone], extraVarsMonotone)
		extraVarsNonmonotone <- data.frame(names(inputDataList)[idxNonmonotone], extraVarsNonmonotone)
		colnames(estAICNonmonotone)[1] <- colnames(inputDataset)[colID]
		colnames(estEDpNonmonotone)[1] <- colnames(inputDataset)[colID]
		colnames(estAICMonotone)[1] <- colnames(inputDataset)[colID]
		colnames(estEDpMonotone)[1] <- colnames(inputDataset)[colID]
		colnames(extraVarsMonotone) <- c(colnames(inputDataset)[colID], paste(all.vars(as.formula(addCovars)), "- MinAIC (est.)"),
				paste(all.vars(as.formula(addCovars)), "- MinAIC (StdErr.)"),
				paste(all.vars(as.formula(addCovars)), "- modelAveraging (est.)"),
				paste(all.vars(as.formula(addCovars)), "- modelAveraging (StdErr.)"))
		colnames(extraVarsNonmonotone) <- c(colnames(inputDataset)[colID], paste(all.vars(as.formula(addCovars)), "- MinAIC (est.)"),
				paste(all.vars(as.formula(addCovars)), "- MinAIC (StdErr.)"),
				paste(all.vars(as.formula(addCovars)), "- modelAveraging (est.)"),
				paste(all.vars(as.formula(addCovars)), "- modelAveraging (StdErr.)"))
		if (!addCovarsVar){
			extraVarsMonotone <- extraVarsMonotone[,-c(4, 5, 8, 9)]
			extraVarsNonmonotone <- extraVarsNonmonotone[,-c(4, 5, 8, 9)]
		}
	}else{
		for (iMono in seq_along(idxMonotone)){
			fittedModelsMonotone[[iMono]] <- resultsMonotone[[iMono]]$fittedModels
			estAICMonotone[iMono,] <- resultsMonotone[[iMono]]$estAIC
			estEDpMonotone[iMono,] <- resultsMonotone[[iMono]]$estEDp
		}
		colnames(estAICMonotone) <- funcListMonotone
		colnames(estEDpMonotone) <- c(funcListMonotone, "minAIC", "modelAveragingAIC")
		for (iNonmono in seq_along(idxNonmonotone)){
			fittedModelsNonmonotone[[iNonmono]] <- resultsNonmonotone[[iNonmono]]$fittedModels
			estAICNonmonotone[iNonmono,] <- resultsNonmonotone[[iNonmono]]$estAIC
			estEDpNonmonotone[iNonmono,] <- resultsNonmonotone[[iNonmono]]$estEDp
		}
		colnames(estAICNonmonotone) <- funcListNonmonotone
		colnames(estEDpNonmonotone) <- c(funcListNonmonotone, "minAIC", "modelAveragingAIC")
		## adding ID's
		fittedModelsFinal<- as.list(rep(NA, length(inputDataList)))
		fittedModelsFinal[idxMonotone] <- fittedModelsMonotone
		fittedModelsFinal[idxNonmonotone] <- fittedModelsNonmonotone
		estAICNonmonotone <- data.frame(names(inputDataList)[idxNonmonotone],estAICNonmonotone)
		estEDpNonmonotone <- data.frame(names(inputDataList)[idxNonmonotone],estEDpNonmonotone)
		estAICMonotone <- data.frame(names(inputDataList)[idxMonotone],estAICMonotone)
		estEDpMonotone <- data.frame(names(inputDataList)[idxMonotone],estEDpMonotone)
		colnames(estAICNonmonotone)[1] <- colnames(inputDataset)[colID]
		colnames(estEDpNonmonotone)[1] <- colnames(inputDataset)[colID]
		colnames(estAICMonotone)[1] <- colnames(inputDataset)[colID]
		colnames(estEDpMonotone)[1] <- colnames(inputDataset)[colID]
	}
	
	

	toReturn <- list(fittedModels = fittedModelsFinal, estAICNonmonotone = estAICNonmonotone, 
			estEDpNonmonotone = estEDpNonmonotone, estAICMonotone = estAICMonotone, 
			estEDpMonotone = estEDpMonotone, extraCovarsMonotone = extraVarsMonotone, 
			extraCovarsNonmonotone = extraVarsNonmonotone)
	class(toReturn) <- "fittedDRM"
	return(toReturn)
}

