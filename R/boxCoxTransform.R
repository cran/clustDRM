# TODO: Add comment
# 
# Author: Vahid Nassiri
###############################################################################
# #' @importFrom caret BoxCoxTrans



#' computes Box-Cox tranform of the response
#' @param x 
#' @return vector the transfomed variable
#' @importFrom caret BoxCoxTrans
#' @author Vahid Nassiri
#' @noRd
boxCoxTransform <- function(x){
	lambda <- BoxCoxTrans(x)$lambda  
	if(lambda==0){
		y=log(x)
	}else{
		y=((x^lambda)-1)/lambda
	}
	return(y)
}
