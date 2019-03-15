# TODO: Add comment
# 
# Author: Vahid Nassiri
###############################################################################


#' performs LRT 
#' @param clusteredData 
#' @param doseLevels 
#' @param nPermute 
#' @param numReplications 
#' @param colID 
#' @param adjustMethod 
#' @param useSeed 
#' @param nCores 
#' @return 
#' 
#' @import doParallel
#' @import parallel
#' @import foreach
#' @import ORIClust
#' @author Vahid Nassiri
#' @noRd
patternPermuteLRT <- function(clusteredData, doseLevels, nPermute, numReplications, colID,
		adjustMethod = c("BH", "holm", "hochberg", "hommel", "bonferroni", "BY", "fdr", "none"), useSeed, nCores){
	unadjustedPvalues <- rep(NA, nrow(clusteredData))
	
	i <- iLRT <- iMCT <- NULL
	if (nCores > parallel::detectCores()){
		stop(paste("Your system only has ", parallel::detectCores(), "you cannot specify an nCores larger than that."))
	}else if(nCores == 1){
		unadjustedPvalues[which(clusteredData$clusteringResult == "decreasing")] <- apply(clusteredData[which(clusteredData$clusteringResult == "decreasing"),-c(1, colID+1)], 1, patternPermuteLRTMonotoneCompletePerSubject, doseLevels = doseLevels, 
				profilePattern = "decreasing", nPermute = nPermute, numReplications = numReplications, useSeed = useSeed)
		
		unadjustedPvalues[which(clusteredData$clusteringResult == "increasing")] <- apply(clusteredData[which(clusteredData$clusteringResult == "increasing"),-c(1, colID+1)], 1, patternPermuteLRTMonotoneCompletePerSubject, doseLevels = doseLevels, 
				profilePattern = "increasing", nPermute = nPermute, numReplications = numReplications, useSeed = useSeed)
		
		unadjustedPvalues[which(clusteredData$clusteringResult == "complete")] <- apply(clusteredData[which(clusteredData$clusteringResult == "complete"),-c(1, colID+1)], 1, patternPermuteLRTMonotoneCompletePerSubject, doseLevels = doseLevels, 
				profilePattern = "complete.profile", nPermute = nPermute, numReplications = numReplications, useSeed = useSeed)
		idxDownOrUp <- which(clusteredData$clusteringResult!= "flat" & clusteredData$clusteringResult!= "decreasing"
						& clusteredData$clusteringResult!= "increasing" & clusteredData$clusteringResult!= "complete")
		
		linputData <- lapply(seq_len(nrow(clusteredData[idxDownOrUp,-c(1, colID+1)])), function(i) clusteredData[idxDownOrUp,-c(1, colID+1)][i,])
		lclusteringResult <- split(as.matrix(clusteredData$clusteringResult[idxDownOrUp]), row(as.matrix(clusteredData$clusteringResult[idxDownOrUp])))
		ldoseLevels <- rep(list(doseLevels), length(linputData))
		lnPermute <- rep(list(nPermute), length(linputData))
		lnumReplications <- rep(list(numReplications), length(linputData))
		luseSeed <- rep(list(useSeed), length(linputData))
		unadjustedPvalues[idxDownOrUp] <- mapply(patternPermuteLRTUpOrDownPerSubject, linputData, 
				ldoseLevels, lnPermute, lnumReplications, lclusteringResult, luseSeed)
		
	}else{
		if (Sys.info()['sysname'] == "Windows"){
			cl <- makeCluster(nCores, type="PSOCK")  
		}else{
			cl <- makeCluster(nCores, type="FORK")  
		}
		registerDoParallel(cl)  
		clusterSetRNGStream(cl = cl, iseed = useSeed)
		unadjustedPvalues[which(clusteredData$clusteringResult == "decreasing")] <- unlist(foreach(iLRT=1:length(which(clusteredData$clusteringResult == "decreasing"))) %dopar% 
						patternPermuteLRTMonotoneCompletePerSubject(clusteredData[which(clusteredData$clusteringResult == "decreasing")[iLRT],-c(1, colID+1)], doseLevels = doseLevels, 
								profilePattern = "decreasing", nPermute = nPermute, numReplications = numReplications, useSeed = useSeed))
		
		unadjustedPvalues[which(clusteredData$clusteringResult == "increasing")] <- unlist(foreach(iLRT=1:length(which(clusteredData$clusteringResult == "increasing"))) %dopar% 
						patternPermuteLRTMonotoneCompletePerSubject(clusteredData[which(clusteredData$clusteringResult == "increasing")[iLRT],-c(1, colID+1)], doseLevels = doseLevels, 
								profilePattern = "increasing", nPermute = nPermute, numReplications = numReplications, useSeed = useSeed))
		
		unadjustedPvalues[which(clusteredData$clusteringResult == "complete")] <- unlist(foreach(iLRT=1:length(which(clusteredData$clusteringResult == "complete"))) %dopar% 
						patternPermuteLRTMonotoneCompletePerSubject(clusteredData[which(clusteredData$clusteringResult == "complete")[iLRT],-c(1, colID+1)], doseLevels = doseLevels, 
								profilePattern = "complete.profile", nPermute = nPermute, numReplications = numReplications, useSeed = useSeed))
		
		idxDownOrUp <- which(clusteredData$clusteringResult!= "flat" & clusteredData$clusteringResult!= "decreasing"
						& clusteredData$clusteringResult!= "increasing" & clusteredData$clusteringResult!= "complete")
		
		#linputData <- lapply(seq_len(nrow(clusteredData[idxDownOrUp,-c(1, colID+1)])), function(i) clusteredData[idxDownOrUp,-c(1, colID+1)][i,])
		lclusteringResult <- split(as.matrix(clusteredData$clusteringResult[idxDownOrUp]), row(as.matrix(clusteredData$clusteringResult[idxDownOrUp])))
		#ldoseLevels <- rep(list(doseLevels), length(linputData))
		#lnPermute <- rep(list(nPermute), length(linputData))
		#lnumReplications <- rep(list(numReplications), length(linputData))
		#luseSeed <- rep(list(useSeed), length(linputData))
		unadjustedPvalues[idxDownOrUp] <- unlist(foreach(iMCT = 1:length(idxDownOrUp)) %dopar% patternPermuteLRTUpOrDownPerSubject(clusteredData[idxDownOrUp,-c(1, colID+1)][iMCT,],
						doseLevels, nPermute, numReplications, lclusteringResult[[iMCT]], useSeed))
		
		
		stopCluster(cl)
	}
	#	unadjustedPvalues[which(clusteredData$clusteringResult == "decreasing")] <- apply(clusteredData[which(clusteredData$clusteringResult == "decreasing"),-c(1, colID+1)], 1, patternPermuteLRTMonotoneCompletePerSubject, doseLevels = doseLevels, 
#			profilePattern = "decreasing", nPermute = nPermute, numReplications = numReplications, useSeed = useSeed)
#	
#	unadjustedPvalues[which(clusteredData$clusteringResult == "increasing")] <- apply(clusteredData[which(clusteredData$clusteringResult == "increasing"),-c(1, colID+1)], 1, patternPermuteLRTMonotoneCompletePerSubject, doseLevels = doseLevels, 
#			profilePattern = "increasing", nPermute = nPermute, numReplications = numReplications, useSeed = useSeed)
#	
#	unadjustedPvalues[which(clusteredData$clusteringResult == "complete")] <- apply(clusteredData[which(clusteredData$clusteringResult == "complete"),-c(1, colID+1)], 1, patternPermuteLRTMonotoneCompletePerSubject, doseLevels = doseLevels, 
#			profilePattern = "complete.profile", nPermute = nPermute, numReplications = numReplications, useSeed = useSeed)
#	idxDownOrUp <- which(clusteredData$clusteringResult!= "flat" & clusteredData$clusteringResult!= "decreasing"
#	& clusteredData$clusteringResult!= "increasing" & clusteredData$clusteringResult!= "complete")
#	
#
#
#	linputData <- lapply(seq_len(nrow(clusteredData[idxDownOrUp,-c(1, colID+1)])), function(i) clusteredData[idxDownOrUp,-c(1, colID+1)][i,])
#	lclusteringResult <- split(as.matrix(clusteredData$clusteringResult[idxDownOrUp]), row(as.matrix(clusteredData$clusteringResult[idxDownOrUp])))
#	ldoseLevels <- rep(list(doseLevels), length(linputData))
#	lnPermute <- rep(list(nPermute), length(linputData))
#	lnumReplications <- rep(list(numReplications), length(linputData))
#	luseSeed <- rep(list(useSeed), length(linputData))
#	unadjustedPvalues[idxDownOrUp] <- mapply(patternPermuteLRTUpOrDownPerSubject, linputData, 
#			ldoseLevels, lnPermute, lnumReplications, lclusteringResult, luseSeed)
#	
	# computing adjusted p-values
	adjustPvalues <- unadjustedPvalues
	adjustPvalues[!is.na(adjustPvalues)] <- p.adjust(unadjustedPvalues[!is.na(adjustPvalues)], method = adjustMethod)
	return(list(unadjustedPvalues = unadjustedPvalues, adjustPvalues = adjustPvalues))
	
	}
	
