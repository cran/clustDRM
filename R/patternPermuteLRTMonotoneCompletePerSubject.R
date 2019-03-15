# TODO: Add comment
# 
# Author: Vahid Nassiri
###############################################################################




#' performs LRT for monotone and complete patterns
#' @param inputDataSubject 
#' @param doseLevels 
#' @param profilePattern 
#' @param nPermute 
#' @param numReplications 
#' @param useSeed 
#' @return 
#' 
#' @import ORIClust
#' @author Vahid Nassiri
#' @noRd
patternPermuteLRTMonotoneCompletePerSubject <- function(inputDataSubject, doseLevels,
		profilePattern = c("decreasing", "increasing", "complete.profile"), 
		nPermute, numReplications, useSeed){
	repeatedDose <- rep(doseLevels, numReplications)
	# making all indexes
	if (!is.null(useSeed)){
		set.seed(useSeed)
	}
	
	indexAll <- apply(matrix(rep(1:length(inputDataSubject), nPermute),
					length(inputDataSubject), nPermute), 2, sample, size = length(inputDataSubject))
	ObservedStat <- computeLRTestStatMonotoneComplete(inputDataSubject, doseLevels, 
			profilePattern, nPermute, numReplications)
	PermutedStat <- apply(matrix(inputDataSubject[indexAll],nrow(indexAll), ncol(indexAll)), 2, computeLRTestStatMonotoneComplete, doseLevels, 
			profilePattern, nPermute, numReplications)
	estimatedPvalue <- (1 + sum(abs(PermutedStat)>abs(ObservedStat)))/(nPermute+1)
	return(estimatedPvalue)
}
