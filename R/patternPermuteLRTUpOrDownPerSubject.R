# TODO: Add comment
# 
# Author: Vahid Nassiri
###############################################################################



#' performs LRT for non-flat profiles per subject
#' @param inputDataSubject 
#' @param doseLevels 
#' @param nPermute 
#' @param numReplications 
#' @param clusteringResult 
#' @param useSeed 
#' @return 
#' 
#' @import ORIClust
#' @importFrom readr parse_number
#' @author Vahid Nassiri
#' @noRd
patternPermuteLRTUpOrDownPerSubject <- function(inputDataSubject, doseLevels,
		nPermute, numReplications, clusteringResult, 
		useSeed){
	repeatedDose <- rep(doseLevels, numReplications)
	# obtain the first word in the clusteringResult to determine the pattern.
	firstWordClusteringResult <- gsub("([A-Za-z]+).*", "\\1", clusteringResult)
	
	if(firstWordClusteringResult == "down"){
		profilePattern = c("down.up")
	} else if(firstWordClusteringResult == "up"){
		profilePattern = c("up.down")
	} else {
		stop("Please give a clusteringResult beginning with words 'up' or 'down'.")
	}
	
	# obtain the numric part of clusteringResult which is then the maximum (or minimum).
	profileMaxOrMin <- parse_number(clusteringResult)
	# making all indexes
	if (!is.null(useSeed)){
		set.seed(useSeed)
	}
	indexAll <- apply(matrix(rep(1:length(inputDataSubject), nPermute),
					length(inputDataSubject), nPermute), 2, sample, size = length(inputDataSubject))
	ObservedStat <- computeLRTestStatUpOrDown(inputDataSubject, doseLevels, 
			profilePattern, nPermute, numReplications, profileMaxOrMin)
	PermutedStat <- apply(matrix(inputDataSubject[indexAll],nrow(indexAll), ncol(indexAll)), 2, computeLRTestStatUpOrDown, doseLevels, 
			profilePattern, nPermute, numReplications, profileMaxOrMin)
	estimatedPvalue <- (1 + sum(abs(PermutedStat)>abs(ObservedStat)))/(nPermute+1)
	# computing adjusted p-values
	return(estimatedPvalue)
}
