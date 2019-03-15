# TODO: Add comment
# 
# Author: Vahid Nassiri
###############################################################################


#' performs MCT for non-flat profiles per subject
#' @param inputDataSubject 
#' @param repeatedDose 
#' @return 
#' 
#' @import multcomp
#' @import doParallel
#' @author Vahid Nassiri
#' @noRd

nonFlatMCTperSubject <- function(inputDataSubject, repeatedDose){
	anovaOnSubject <- aov(unlist(inputDataSubject) ~ repeatedDose)
	glhtModel <- summary(multcomp::glht(anovaOnSubject,linfct = multcomp::mcp(repeatedDose = "UmbrellaWilliams")))
	MCTpvalues <- as.vector(glhtModel$test$pvalues)
	MCTtestStat <- glhtModel$test$tstat
	toReturn <- c(which.max(MCTtestStat), max(MCTtestStat),min(MCTpvalues))	
	names(toReturn) <- c("selectedContrast", "testStat", "pValue")
	return(toReturn)
}

