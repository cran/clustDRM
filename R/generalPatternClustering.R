# TODO: Add comment
# 
# Author: Vahid Nassiri
###############################################################################


#' Clustering dose-response curves based on their pattern
#' 
#' function to cluster dose-response curves based on their pattern.
#' @param inputData data matrix which should incluide ID's of the subjects, as well as the measurements (gene expressions, etc.) 
#' for all replications of different as columns. 
#' @param colsData vector indicating the idex of columns in the inputData which correspond to the measurement for different replications 
#' of different doses.
#' @param colID scalar indicating the index of column corresponding to data ID.
#' @param doseLevels vector with dose levels.
#' @param numReplications vector wit hthe same length as doseLevels with number of replications for each dose. 
#' @param na.rm logical variable indicatign whether missing values should be removed (TRUE) or not (FALSE, default)
#' @param imputationMethod signle string taking calues from "mean" (default), and "median", which indicates how the missing values should be
#' treated. "mean" would replace them with the mean of the observed ones, and "median" will use median of them for imputation.
#' @param ORICC signle string taking value "two", "one", and "both", indicating which ORICC procedure should be used. "one" refers to 
#' one-stage ORICC only, "two" (default) refers to two-stage ORICC only, and "both" will perform both of them.
#' @param transform single string indicating what kind of transform should be applied on the response data. 
#' It takes "none" (no transform, dafault), "log" (natural log), "sRoot (square root), and "qRoot" (cubic root), and 
#' "boxcox" (Box-Cox transformation).
#' @param plotFormat plotFormat string gets two values "eps" (default), and "jpg" indicating the format of the ouput plot.
#' @param LRT logical indicating whether a permutation-based likelihood ratio test should be applied (TRUE) on the subjects which 
#' their trend is identified as non-flast by ORICC1 or not (FALSE).
#' @param MCT logical indicating whether a multiple comparison test (with "UmbrellaWilliams" constrast matrix) 
#' should be applied (TRUE) on the subjects which their trend is identified as non-flast by ORICC1 or not (FALSE).
#' @param adjustMethod The method for multiplicity adjustment for p-values. The possible values for this argument are 
#' "BH", "holm", "hochberg", "hommel", "bonferroni", "BY", "fdr", "none" with "BH" (Benjamini-Hocheberg) as default. 
#' @param nPermute scalar indicating number of permutations in LRT.
#' @param useSeed scalar, indicating the seed should be used to generate LRT permutations. The default is NULL.
#' @param nCores nCores scalar, indicating the number of cores should be used to perform LRT and MCT tests. Default is 1 which means sequantial 
#' computation (no prallel computation).
#' @param theLeastNumberOfMethods scalar taking values from 1, 2, 3, and 4, indicating how many methods should 
#' approve a non-flat trend that it can be selected. Its value depends on how many tests are asked to be done,
#' for the maximum happens when ORICC = "both" and both LRT and MCT are TRUE. For example, when this argument sets to 2 and ORICC = "two", 
#' LRT = TRUE, and MCT = TRUE, it means if two-stage ORICC identifies a non-flat pattern and at least one of the 
#' LRT and MCT also accepts (at the level of alpha), then that comound is selected as one with a non-flat pattern. 
#' Note that the comparison with alpha is done for adjusted p-values.
#' @param alpha the significance level to compare the adjusted p-value with it.
#' @details This function first use ORIIC1 or ORICC2 (or both) to identify the pattern of the dose-response cruve for each subject. Once 
#' the pattern is identified, for non-flat ones, a permutation-based likelihood ratio test (for exactly the identified pattern, 
#' if LRT = TRUE), and a multiple comparisons test (to test H0: flat vs. H1: non-flat, if MCT = TRUE) will be performed to further filter 
#' the flat patterns. 
#' 
#' @return a list of the following objects: 
#' 
#' selectedSubjects: a data frame indicating the ID's of the selected subjects in the first columns and the 
#' identified trend in the second column.
#' 
#' clusteringORICC1Results and/or clusteringORICC2Results: a list with four elements providing the raw data as 
#' the outcome of the ORICC procedure (rawDataORICC1 and/or rawDataORICC2), the pattern identified by the ORICC 
#' procedure (clusteringResultsORICC1 and/or clusteringResultsORICC2), results of LRT (resultsLRT) and results of 
#' MCT (resultsMCT). Both of them provide the adjusted and unadjusted p-values, but for MCT the selected contrast
#' will be provided as well.
#' 
#' @examples
#' ## gnerating data
#' set.seed(11)
#' doses2Use <-  c(0, 5, 20)
#' numRep2Use <- c(3, 3, 3)
#' generatedData <- cbind(rep(1,sum(numRep2Use)),
#' MCPMod::genDFdata("logistic",c(5, 3, 10, 0.05), doses2Use, 
#' numRep2Use, 1), 
#' 		matrix(rnorm(1*sum(numRep2Use)), sum(numRep2Use), 1))
#' colnames(generatedData) <- c("ID", "dose", "response", "x1")
#' for (iGen in 2:15){
#' 	genData0 <- cbind(rep(iGen,sum(numRep2Use)),
#' MCPMod::genDFdata("logistic",c(5, 3, 10, 0.05), doses2Use, 
#' numRep2Use, 1), 
#' 			matrix(rnorm(1*sum(numRep2Use)), sum(numRep2Use), 1))
#' 	colnames(genData0) <- c("ID", "dose", "response", "x1")
#' 	generatedData <- rbind(generatedData, genData0)
#' }
#' ## transforming it for clustering
#' toInput <- inputDataMaker(2, 3, 1, generatedData)
#' ## general pattern clustering
#' generalPatternClust <- generalPatternClustering(inputData = toInput$inputData, 
#' colsData = toInput$colsData ,colID = toInput$colID , 
#' 		doseLevels = toInput$doseLevels, numReplications = toInput$numReplicates, 
#' na.rm = FALSE, imputationMethod = "mean",
#' 		ORICC = "two", transform = "none",plotFormat = "eps", 
#' LRT = TRUE, MCT = TRUE,
#' 		adjustMethod = "BH",
#' 		nPermute = 100, useSeed = NULL, 
#' theLeastNumberOfMethods = 2, alpha = 0.05, nCores = 1)
#' @import doParallel
#' @import ORIClust
#' @import multcomp
#' @import parallel
#' @importFrom parallel detectCores
#' @importFrom MCPMod genDFdata
#' @seealso \href{https://www.rdocumentation.org/packages/ORIClust/versions/1.0-1/topics/ORICC1}{ORIClust}
#' \href{https://www.rdocumentation.org/packages/ORIClust/versions/1.0-1/topics/ORICC2}{ORIClust}
#' @author Vahid Nassiri, and Yimer Wasihun.
#' @export

generalPatternClustering <- function(inputData, colsData ,colID , doseLevels, numReplications, na.rm = FALSE, imputationMethod = c("mean", "median"),
		ORICC = c("two", "one", "both"), transform = c("none", "log","sRoot", "qRoot", "boxcox"), plotFormat = c("eps", "jpg"), LRT = TRUE, MCT = FALSE,
		adjustMethod = c("BH", "holm", "hochberg", "hommel", "bonferroni", "BY", "fdr", "none"),
		nPermute = 1000, useSeed = NULL, theLeastNumberOfMethods = c(1, 2, 3, 4), alpha = 0.05, nCores = 1){
	
	## cheking arguments
	if (colID<1 | colID>ncol(inputData)){
		stop(paste("colID should be within 1 and ", ncol(inputData)))
	}
	if (min(colsData) < 1 | max(colsData) > ncol(inputData)){
		stop("Enter colsData within 1 and ", ncol(inputData))
	}
	if (length(doseLevels) != length(numReplications)){
		stop("doseLevels should be of the same length as numReplications.")
	}
	imputationMethod <- match.arg(imputationMethod,  c("mean", "median"))
	ORICC <- match.arg(ORICC,  c("two", "one", "both"))
	transform <- match.arg(transform, c("none", "log","sRoot", "qRoot", "boxcox"))
	plotFormat <- match.arg(plotFormat,  c("eps", "jpg"))
	adjustMethod <- match.arg(adjustMethod,  c("BH", "holm", "hochberg", "hommel", "bonferroni", "BY", "fdr", "none"))
	if (nPermute < 1){
		stop("nPermute should be larger than 0.")
	}
	if (nCores > parallel::detectCores()){
		stop(paste("Your system only has ", parallel::detectCores(), "you cannot specify nCores larger than that."))
	}
	if (ORICC == "both" & theLeastNumberOfMethods > (2 + LRT + MCT)){
		stop("theLeastNumberOfMethods cannot be larger than ", 2 + LRT + MCT)
	}
	if (ORICC != "both" & theLeastNumberOfMethods > (1 + LRT + MCT)){
		stop("theLeastNumberOfMethods cannot be larger than ", 1 + LRT + MCT)
	}
	if (alpha <0 | alpha >1){
		stop("alpha should be in (0,1).")
	}
	## keeping original input data
	#origInputData <- inputData
	## As we are constructing that outseveles, so we know it's in column 3
	## finding subject with missing values
	idxNA <- which(is.na(inputData), arr.ind=TRUE)
	idxNARows <- unique(idxNA[,1])
	## either remove subjects with missing values or replace them with
	## the mean or median of the avilable values.
	clusteringORICC1Results <- NULL
	clusteringORICC2Results <- NULL
	if (na.rm){
		inputData <- inputData[-idxNARows,]
	}else{
		missingSubjects <- inputData[idxNARows,-colID]
		imputedValues <- apply(missingSubjects, 1, get(imputationMethod), na.rm = TRUE)
		for (iNA in 1:length(idxNARows)){
			inputData[idxNARows[iNA], idxNA[idxNA[,1]== idxNARows[iNA],2]] <- imputedValues[iNA]
		}
	}
	## transforming data if needed
	if (transform == "none"){
		inputData[,colsData] <- inputData[,colsData]
	}else if(transform == "log"){
		if (any(unlist(inputData[,colsData]) == 0 )){
			inputData[,colsData] <- log(inputData[,colsData] + 1)
			warning("Due to zero responses, log(x+1) is used instead of log.")
		}
		if(any(unlist(inputData[,colsData]) <0)){
			stop("log transofrm is not possible with negative responses.")
		}
		if(all(unlist(inputData[,colsData])>0)){
			inputData[,colsData] <- log(inputData[,colsData])
		}
	}else if (transform == "sRoot"){
		if (any(unlist(inputData[,colsData])<0)){
			stop("square root is not possible with negative responses.")
		}else{
			inputData[,colsData] <- sqrt(inputData[,colsData])
		}
	}else if (transform == "qRoot"){
		inputData[,colsData] <- sign(inputData[,colsData]) * (abs(inputData[,colsData])^(1/3))
	}else if(transform  == "boxcox"){
		inputData[,colsData] <- apply(inputData[,colsData], 2, boxCoxTransform)
	}

	## performing clustering
	if (ORICC == "two"){
		## clustering
		clusteringORICC2Results <- clusteringORICC2 (inputData, colsData, colID, doseLevels, numReplications, 
				adjustMethod, LRT, MCT, nPermute, useSeed, transform = 0, plotFormat, nCores)
		
		clusteringSummary <- apply(cbind(clusteringORICC2Results$clusteringResultsORICC2 != "flat", clusteringORICC2Results$resultsLRT$`adjusted p-value` < alpha, 
						clusteringORICC2Results$resultsMCT$adjustedPvalues < alpha),1,sum, na.rm = TRUE)
		idxNonflat <- which(clusteringSummary >= theLeastNumberOfMethods)
		selectedSubjects <- data.frame(inputData[idxNonflat,colID],
				clusteringORICC2Results$clusteringResultsORICC2[idxNonflat])
		colnames(selectedSubjects) <- c(names(inputData)[colID], "pattern")
	}else if(ORICC == "one"){
		## clustering
		clusteringORICC1Results <- clusteringORICC1 (inputData, colsData, colID, doseLevels, 
				numReplications, adjustMethod, LRT, MCT, nPermute, useSeed, transform = 0, plotFormat, nCores)
		
		clusteringSummary <- apply(cbind(clusteringORICC1Results$clusteringResultsORICC1 != "flat", clusteringORICC1Results$resultsLRT$`adjusted p-value` < alpha, 
						clusteringORICC1Results$resultsMCT$adjustedPvalues < alpha),1,sum, na.rm = TRUE)
		idxNonflat <- which(clusteringSummary >= theLeastNumberOfMethods)
		selectedSubjects <- data.frame(inputData[idxNonflat,colID],
				clusteringORICC1Results$clusteringResultsORICC1[idxNonflat])
		colnames(selectedSubjects) <- c(names(inputData)[colID], "pattern")
				
	}	else if(ORICC == "both"){
		## clustering
		clusteringORICC2Results <- clusteringORICC2 (inputData, colsData, colID, doseLevels, 
				numReplications, adjustMethod, LRT, MCT, nPermute, useSeed, transform = 0, plotFormat, nCores)
		clusteringORICC1Results <- clusteringORICC1 (inputData, colsData, colID, doseLevels, 
				numReplications, adjustMethod, LRT, MCT, nPermute, useSeed, transform = 0, plotFormat, nCores)
		
		#class(clusteringORICC1Results) <- "nonmonoClust"
		#class(clusteringORICC2Results) <- "nonmonoClust"
	clusteringSummary <- apply(cbind(clusteringORICC1Results$clusteringResultsORICC1 != "flat", clusteringORICC1Results$resultsLRT$`adjusted p-value` < alpha, 
				clusteringORICC1Results$resultsMCT$adjustedPvalues < alpha,
				clusteringORICC2Results$clusteringResultsORICC2 != "flat", clusteringORICC2Results$resultsLRT$`adjusted p-value` < alpha, 
				clusteringORICC2Results$resultsMCT$adjustedPvalues < alpha),1,sum, na.rm = TRUE)

	idxNonflat <- which(clusteringSummary >= theLeastNumberOfMethods)
	selectedSubjects <- data.frame(inputData[idxNonflat,colID],
			clusteringORICC1Results$clusteringResultsORICC1[idxNonflat])
	colnames(selectedSubjects) <- c(names(inputData)[colID], "pattern")


	}
	toReturn <- list(selectedSubjects = selectedSubjects, clusteringORICC1Results = clusteringORICC1Results, clusteringORICC2Results = clusteringORICC2Results)
	return(toReturn)	
}

