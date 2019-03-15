# TODO: Add comment
# 
# Author: Vahid Nassiri
###############################################################################


#' performs MCT for non-flat profiles
#' @param clusteredData 
#' @param doseLevels 
#' @param numReplications 
#' @param colID 
#' @param adjustMethod 
#' @param nCores 
#' @return 
#' 
#' @import doParallel
#' @import parallel
#' @import foreach
#' @author Vahid Nassiri
#' @noRd
nonFlatMCT <- function(clusteredData, doseLevels, numReplications, colID,
		adjustMethod = c("BH", "holm", "hochberg", "hommel", "bonferroni", "BY", "fdr", "none"), nCores){
	i <- iLRT <- iMCT <- NULL
	
	compoundData <- clusteredData[, -c(1, colID+1)]
	repeatedDose <- as.factor(rep(doseLevels, numReplications))
	idxNonFlat <- which(clusteredData$clusteringResult != "flat")
	
	
	if (nCores > parallel::detectCores()){
		stop(paste("Your system only has ", parallel::detectCores(), "you cannot specify an nCores larger than that."))
	}else if(nCores == 1){
		resultsMCT <- apply(compoundData[idxNonFlat,], 1, nonFlatMCTperSubject, repeatedDose)
	}else{
		if (Sys.info()['sysname'] == "Windows"){
			cl <- makeCluster(nCores, type="PSOCK")  
		}else{
			cl <- makeCluster(nCores, type="FORK")  
		}
		registerDoParallel(cl)  
		resultsMCTList <- foreach(i=1:length(idxNonFlat)) %dopar% nonFlatMCTperSubject(compoundData[idxNonFlat[i],], repeatedDose)
		stopCluster(cl)
		resultsMCT <- matrix(unlist(resultsMCTList), 3, length(resultsMCTList))
	}
	
	unadjustedPvalues <- rep(NA, nrow(clusteredData))
	testStat <- rep(NA, nrow(clusteredData))
	selectedContrast <- rep(NA, nrow(clusteredData))
	unadjustedPvalues[idxNonFlat] <- resultsMCT[3,]
	testStat[idxNonFlat] <- resultsMCT[2,]
	selectedContrast[idxNonFlat] <- resultsMCT[1,]
	adjustedPvalues <- unadjustedPvalues
	adjustedPvalues[idxNonFlat] <- p.adjust(unadjustedPvalues[idxNonFlat], method = adjustMethod)
	return(list(testStat = testStat, unadjustedPvalues = unadjustedPvalues, adjustedPvalues = adjustedPvalues, 
					selectedContrast = selectedContrast))	
}
