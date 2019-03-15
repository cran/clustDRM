# TODO: Add comment
# This function is obtained from calcED function in DoseFinding package 
# (https://cran.r-project.org/web/packages/DoseFinding/)
# Author: Vahid Nassiri
###############################################################################


#' computed effective dose for a given set of parameters
#' @param model 
#' @param pars 
#' @param p 
#' @param maxD 
#' @param EDtype 
#' @param doses 
#' @param off 
#' @param scal 
#' @param nodes 
#' @return 
#' 
#' @author Vahid Nassiri (taken from DoseFinding package)
#' @noRd

localCalcED <- function(model, pars, p, maxD, EDtype = c("continuous", "discrete"),
		doses, off, scal, nodes){
	## calculate the smallest dose x for which
	## f(x) > f(0) + p*(f(xmax)-f(0))
	## e.g. the EDp within the observed dose-range
	EDtype <- match.arg(EDtype)
	if(model == "betaMod" & missing(scal))
		stop("Need \"scal\" parameter for betaMod model")
	if(model == "linlog" & missing(off))
		stop("Need \"off\" parameter for linlog model")    
	if(model == "linInt"){
		if(missing(nodes))
			stop("Need \"nodes\" parameter for linlog model")
		if(length(nodes) != length(pars))
			stop("nodes and pars of incompatible length")
	}
	
	if(EDtype == "continuous"){ ## calculate target dose analytically
		cf <- pars
		if(cf[2] == 0){
			return(NA)
		}
		if(model == "linear"){
			return(p*maxD)
		}
		if(model == "linlog"){
			return(off*(exp(p*(log(maxD+off)-log(off)))-1))
		}
		if(model == "exponential"){
			return(cf[3]*log(p*exp(maxD/cf[3])-p+1))
		}
		if(model == "emax"){
			return(p*cf[3]*maxD/((1-p)*maxD+cf[3]))
		}
		if(model == "logistic"){
			res1 <- ((p-1)*exp(maxD/cf[4]+cf[3]/cf[4])-exp(2*cf[3]/cf[4])-p*exp(cf[3]/cf[4]))
			res2 <- ((p*exp(cf[3]/cf[4])+1)*exp(maxD/cf[4])+(1-p)*exp(cf[3]/cf[4]))
			return(cf[3]-cf[4]*log(-res1/res2))
		}
		if(model == "sigEmax"){
			out <-  p*cf[3]^cf[4]*maxD^cf[4]/((1-p)*maxD^cf[4]+cf[3]^cf[4])
			return(out^(1/cf[4]))
		}
		if(model == "quadratic"){
			mode <- -pars[2]/(2*pars[3])
			if(mode > maxD | mode < 0) ## maximum outside dose range
				mode <- maxD
			const <- pars[2]*mode+pars[3]*mode^2
			d1 <- -(sqrt(4*pars[3]*const*p+pars[2]^2)+pars[2])/pars[3]/2.0
			d2 <- (sqrt(4*pars[3]*const*p+pars[2]^2)-pars[2])/pars[3]/2.0
			ind <- c(d1, d2) > 0
			if(!any(ind))
				return(NA)
			return(min(c(d1, d2)[ind]))
		}
		if(model == "betaMod"){
			func <- function(x, Emax, delta1, delta2, scal, p, mode){
				p - betaMod(x, 0, 1, delta1, delta2, scal)/betaMod(mode, 0, 1, delta1, delta2, scal)
			}
			mode <- cf[3]/(cf[3]+cf[4])*scal
			out <- uniroot(func, lower=0, upper=mode, delta1=cf[3],
					delta2=cf[4], Emax=cf[2], scal=scal,
					p=p, mode = mode)$root
			return(out)
		}
	}
	
	
	if(EDtype == "discrete"){
		warning("In this version we do not support discrete ED type, please use Dosefinding.")
	}
	
}
