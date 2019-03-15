# TODO: Add comment
# 
# Author: Vahid Nassiri
###############################################################################


#' performs clustering using ORICC2
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
clusteringORICC2 <- function(inputData, colsData, colID, doseLevels = NULL, numReplications, 
		adjustMethod = c("BH", "holm", "hochberg", "hommel", "bonferroni", "BY", "fdr", "none"), 
		LRT, MCT, nPermute = NULL, useSeed,  
		transform = 0,  plotFormat, nCores){
	
	# dose levels are needed if any of the tests needs to be done
	if ((LRT | MCT) & is.null(doseLevels)){
		stop("doseLevels are needed to perform LRT or MCT tests!")
	}
	if (LRT & is.null(nPermute)){
		stop("nPermute is needed to perform LRT!")
	}
	
	#nameRawDataOutputORICC2 = paste(savingLocation, paste("/",outcomeName,"_ORICC2RawData.txt",sep=""), sep="")
	#nameFittedDataOutputORICC2 = paste(savingLocation, paste("/",outcomeName,"_ORICC2MeanData.txt",sep=""), sep="")
	#namePlotOutputORICC2 = paste(savingLocation, paste("/",outcomeName,"_ORICC2Plot.eps",sep=""), sep="")
	clusteringORICC2 <- modifiedORICC2(data = inputData, data.col = colsData, id.col=colID, n.rep=numReplications, 
			n.top = NULL,transform = 0, name.profile="all", cyclical.profile=NULL, onefile=NULL,
			plot.format = plotFormat)
	
	
	
	rawDataORICC2 <- clusteringORICC2$rawDataWithClustering
	# matching outcome of clustering with subjects for ORICC2
	clusteringResultsORICC2 <- rep('flat', length(inputData[,colID]))
	####################3
	###################### NEED CORRECTION
	#colIDORICC2 <- 3
	colIDORICC2 <- which(names(rawDataORICC2) == names(inputData)[colID])
	matchORICC2 <- match(rawDataORICC2[, colIDORICC2], inputData[,colID])
	clusteringResultsORICC2[matchORICC2] <- as.character(rawDataORICC2$clusterNames)
	
	## now doing LRT if asked
	
	clusteredData2 <- data.frame(clusteringResultsORICC2, inputData)
	names(clusteredData2)[1] <- "clusteringResult"
	outputORICC2 <- clusteredData2
	names(outputORICC2)[2] <- "ID"
	outputLRT2 <- NULL
	if (LRT){
		patternPermuteLRT2 <- patternPermuteLRT(clusteredData = clusteredData2, doseLevels, nPermute, numReplications, colID,
				adjustMethod = adjustMethod, useSeed = useSeed, nCores)
		outputLRT2 <- data.frame(clusteringResultsORICC2, inputData[,colID], patternPermuteLRT2$unadjustedPvalues, patternPermuteLRT2$adjustPvalues)
		names(outputLRT2) <- c("clusteringResult","ID", "unadjusted p-value", "adjusted p-value")
		outputLRT <- list(rawDataORICC2 = rawDataORICC2, clusteringResultsORICC2 = clusteringResultsORICC2, resultsLRT = outputLRT2)
	}
	outputMCT2 <- NULL
	if (MCT){
		nonflatMCT2 <- nonFlatMCT(clusteredData2, doseLevels, numReplications, colID,
				adjustMethod, nCores)
		outputMCT2 <- data.frame(clusteringResultsORICC2, inputData[,colID],nonflatMCT2$testStat, nonflatMCT2$unadjustedPvalues,
				nonflatMCT2$adjustedPvalues, nonflatMCT2$selectedContrast)
		names(outputMCT2) <- c("clusteringResults", "ID", "testStat", "unadjustedPvalues", "adjustedPvalues", "selectedContrast")
	}
	
	
	return(list(rawDataORICC2 = rawDataORICC2, clusteringResultsORICC2 = clusteringResultsORICC2, resultsLRT = outputLRT2, 
					resultsMCT = outputMCT2))
}

