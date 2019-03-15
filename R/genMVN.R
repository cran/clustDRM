# TODO: Add comment
# 
# Author: Vahid Nassiri
###############################################################################


#' generates random numbers from multivariate normal
#' @param n 
#' @param mu 
#' @param sigmaMat 
#' @return 
#' 
#' @author Vahid Nassiri
#' @noRd

genMVN <- function(n, mu, sigmaMat){
	## obtain number of variables
	numVar <- length(mu)
	## generate independent data
	indData <- matrix(rnorm(n*numVar),n , numVar)
	## find the Cholesky decomposition of the covariance matrix
	choleskyDecompSigma <- chol(sigmaMat)
	## make the indepedent data dependent
	genData <- (indData%*% choleskyDecompSigma) + mu[col(indData%*% choleskyDecompSigma)]
	return(genData)
}
