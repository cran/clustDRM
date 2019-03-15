# TODO: Add comment
# 
# Author: Vahid Nassiri
###############################################################################



#' computes LRT test statistic for umbrella-shaped profiles
#' @param inputDataSubject 
#' @param doseLevels 
#' @param profilePattern 
#' @param nPermute 
#' @param numReplications 
#' @param profileMaxOrMin 
#' @return 
#' 
#' @author Vahid Nassiri
#' @noRd

computeLRTestStatUpOrDown <- function(inputDataSubject, doseLevels, 
		profilePattern = c("ORIClust::up.down", "ORIClust::down.up"), 
		nPermute, numReplications, profileMaxOrMin){
	## it is only applied to decreasing, increasing, and complete.profile.
	## note that input data is a vector: the data for a compound, gene, etc.
	## compute mean per dose
	inputDataSubject <- as.numeric(inputDataSubject)
	dose <- rep(doseLevels, numReplications)
	meanPerDose <- tapply(X = inputDataSubject, INDEX = as.factor(dose), FUN = mean, na.rm = TRUE)
	profilePatternFunction <- get(profilePattern)
	profileFit <- profilePatternFunction(inputDataSubject, meanPerDose, numReplications, h = profileMaxOrMin)
	## compute the test statistic
	observedRSS0   <- sum((inputDataSubject-mean(inputDataSubject))^2)
	observedRSS1   <- sum((inputDataSubject-rep(profileFit$mu, numReplications))^2)
	observedStatistic <- 1-(observedRSS1/observedRSS0)
	return(observedStatistic)
}

