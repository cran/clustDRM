

#' performs clustering using ORICC1
#' @param inputData 
#' @param colsData 
#' @param colID 
#' @param doseLevels 
#' @param numReplications 
#' @param adjustMethod 
#' @param LRT 
#' @param MCT 
#' @param nPermute 
#' @param useSeed 
#' @param transform 
#' @param plotFormat 
#' @param nCores 
#' @return 
#' 
#' @author Vahid Nassiri
#' @noRd
clusteringORICC1 <- function(inputData, colsData, colID, doseLevels , numReplications, 
		adjustMethod = c("BH", "holm", "hochberg", "hommel", "bonferroni", "BY", "fdr", "none"), 
		LRT, MCT, nPermute = NULL, useSeed,  
		transform,  plotFormat, nCores){
	# dose levels are needed if any of the tests needs to be done
	if ((LRT | MCT) & is.null(doseLevels)){
		stop("doseLevels are needed to perform LRT or MCT tests!")
	}
	if (LRT & is.null(nPermute)){
		stop("nPermute is needed to perform LRT!")
	}
	clusteringORICC1 <- modifiedORICC1(data = inputData, data.col= colsData, id.col=colID, n.rep=numReplications,
			n.top = NULL,transform = transform, name.profile="all", complete.profile=1, onefile=NULL,
			plot.format= plotFormat)
	rawDataORICC1 <- clusteringORICC1$rawDataWithClustering
	# matching outcome of clustering with subjects for ORICC1
	clusteringResultsORICC1 <- rep('flat', length(inputData[,colID]))
	colIDORICC1 <- which(names(rawDataORICC1) == names(inputData)[colID])
	matchORICC1 <- match(rawDataORICC1[, colIDORICC1], inputData[,colID])
	clusteringResultsORICC1[matchORICC1] <- as.character(rawDataORICC1$clusterNames)
	clusteredData1 <- data.frame(clusteringResultsORICC1, inputData)
	names(clusteredData1)[1] <- "clusteringResult"
	outputORICC1 <- clusteredData1
	names(outputORICC1)[2] <- "ID"
	outputLRT1 <- NULL
	if (LRT){
		patternPermuteLRT1 <- patternPermuteLRT(clusteredData1, doseLevels, nPermute, numReplications, colID,
				adjustMethod = adjustMethod, useSeed = useSeed, nCores)
		outputLRT1 <- data.frame(clusteringResultsORICC1, inputData[,colID], patternPermuteLRT1$unadjustedPvalues, patternPermuteLRT1$adjustPvalues)
		names(outputLRT1) <- c("clusteringResult","ID", "unadjusted p-value", "adjusted p-value")
		outputLRT <- list(rawDataORICC1 = rawDataORICC1, clusteringResultsORICC1 = clusteringResultsORICC1, resultsLRT = outputLRT1)
	}
	outputMCT1 <- NULL
	if (MCT){
		nonflatMCT1 <- nonFlatMCT(clusteredData1, doseLevels, numReplications, colID,
				adjustMethod, nCores)
		outputMCT1 <- data.frame(clusteringResultsORICC1, inputData[,colID],nonflatMCT1$testStat, nonflatMCT1$unadjustedPvalues,
				nonflatMCT1$adjustedPvalues, nonflatMCT1$selectedContrast)
		names(outputMCT1) <- c("clusteringResultsORICC1", names(inputData)[colID], "testStat",
				"unadjustedPvalues", "adjustedPvalues", "selectedContrast")
	}
	return(list(rawDataORICC1 = rawDataORICC1, clusteringResultsORICC1 = clusteringResultsORICC1, resultsLRT = outputLRT1, 
					
					resultsMCT = outputMCT1))
	
}
