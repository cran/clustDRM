# TODO: Add comment
# 
# Author: Vahid Nassiri
###############################################################################

#' plot dose-response curves
#' 
#' function to plot dose-response curves with the possibility of adding lines indicating average response per dose levels. Also, 
#' provided a pattern for the dose-response curve, it can estimate the expected mean values per dose level for the given pattern and
#' add them to the plot.
#' 
#' @param inputDataset a data frame containing the input dataset, it should at least include dose, response, and ID
#' @param dose either a single string or a scalar, indicating the name of the dose column or its index.
#' @param response either a single string or a scalar, indicating the name of the response column or its index.
#' @param ID either a single string or a scalar, indicating the name of the ID column or its index.
#' @param subjectID single input as the same type as given ID column with the ID of the subject to plot.
#' @param xlab single string with default "dose", the label on x axis. 
#' @param ylab single string with default "response", the label on y axis. 
#' @param addMean logical variable indicating whether mean values (connecting with lines) should be plotted or not.
#' @param drcPattern single string showing the idetified pattern using clustering algorithms. The default is NULL. In such case, no extra 
#' line will be added to the plot regarding the estimated means via the identified pattern.
#' @return make a plot.
#' 
#' @details with addMean = TRUE, a line will be added to the plot, connecting the averaged response per dose level. 
#' But when a pattern is provided for the dose-response curve via drcPattern, then a line will be added to the data with the means 
#' estimated assuming the identified pattern. If both addMEan = TRUE and drcPattern != NULL, then two lines will be added to the plot. The line in 
#' purplish-colored with cross signs as points is the averaged response value per dose level, and the bluish-colored line with circled cross signs as points
#' represents the estimated mean based on the pattern.
#' 
#' @import ORIClust
#' @importFrom MCPMod genDFdata 
#' @importFrom graphics boxplot lines par plot points stripchart
#' @examples 
#' ## gnerating data, a sample of size 20
#' set.seed(11)
#' doses2Use <-  c(0, 5, 20)
#' numRep2Use <- c(3, 3, 3)
#' generatedData <- cbind(rep(1,sum(numRep2Use)), 
#' MCPMod::genDFdata("logistic",c(5, 3, 10, 0.05), doses2Use, 
#' numRep2Use, 1), 
#' 		matrix(rnorm(1*sum(numRep2Use)), sum(numRep2Use), 1))
#' colnames(generatedData) <- c("ID", "dose", "response", "x1")
#' for (iGen in 2:15){
#' genData0 <- cbind(rep(iGen,sum(numRep2Use)), 
#' MCPMod::genDFdata("logistic",c(5, 3, 10, 0.05), 
#' doses2Use, numRep2Use, 1), 
#' 			matrix(rnorm(1*sum(numRep2Use)), sum(numRep2Use), 1))
#' colnames(genData0) <- c("ID", "dose", "response", "x1")
#' generatedData <- rbind(generatedData, genData0)
#' }
#' ## plotting dose response relation
#' plotDoseResponseData(generatedData, 2, 3, 1, 2)
#' ## transforming it for clustering
#' plotDoseResponseData(generatedData, 2, 3, 1, 2, 
#' addMean = FALSE, 
#' 		drcPattern = "increasing")
#' 
#' @author Vahid Nassiri and Yimer Wasihun.
#' @export

plotDoseResponseData <- function(inputDataset, dose, response, ID, subjectID, xlab = "Dose", ylab = "Response", 
addMean = TRUE, drcPattern = NULL){
	colID <- ID
	doseColumn <- dose 
	responseCol <- response
	if (is.character(ID)){
		colID <- which(names(inputDataset) == ID)
	}
	if (is.character(dose)){
		doseColumn <- which(names(inputDataset) == dose)
	}
	if (is.character(response)){
		responseCol <- which(names(inputDataset) == response)
	}
	idxSubject <- which(inputDataset[,colID] == subjectID)
	doseData <- inputDataset[,doseColumn][idxSubject]
	responseData <- inputDataset[,responseCol][idxSubject]
	boxplot(unlist(responseData)~as.factor(doseData), xlab = xlab, 
			ylab = ylab)
	stripchart(unlist(responseData)~as.factor(doseData), vertical = TRUE, method = "jitter",
			pch = 21, col = "maroon", bg = "bisque",
			add = TRUE)
	if (addMean){
		data2Plot <- data.frame(unlist(responseData), as.factor(doseData))
		names(data2Plot) <- c("response", "dose")
		# meanValues <- plyr::ddply(data2Plot, ~dose, plyr::summarize, mean(response))
		meanPerDose <- rep(NA, length(levels(data2Plot$dose)))
		for (iDose in 1:length(levels(data2Plot$dose))){
			meanPerDose[iDose] <- mean(data2Plot$response[data2Plot$dose == levels(data2Plot$dose)[iDose]])
		}
		meanValues <- data.frame(levels(data2Plot$dose), meanPerDose)
		names(meanValues)[1] <- "dose"
		points(meanValues$dose, meanValues[,2], pch = 4, col = "maroon4", cex = 1.2)
		lines(meanValues$dose, meanValues[,2], col = "orchid3", lwd = 1.5)
	}
	if (!is.null(drcPattern)){
		if (drcPattern == "increasing"){
			estMu <- ORIClust::increasing(responseData, 
					tapply(X = responseData, INDEX = as.factor(doseData), FUN = mean, na.rm = TRUE), table(doseData))$mu
		}else if(drcPattern == "decreasing"){
			estMu <- ORIClust::decreasing(responseData, 
					tapply(X = responseData, INDEX = as.factor(doseData), FUN = mean, na.rm = TRUE), table(doseData))$mu
		}else if(drcPattern == "flat"){
			estMu <- ORIClust::flat.pattern(responseData, 
					tapply(X = responseData, INDEX = as.factor(doseData), FUN = mean, na.rm = TRUE), table(doseData))$mu
		}else if(gsub("([A-Za-z]+).*", "\\1", drcPattern) == "down"){
			estMu <- ORIClust::down.up(responseData, 
					tapply(X = responseData, INDEX = as.factor(doseData), FUN = mean, na.rm = TRUE), 
					table(doseData), h = parse_number(drcPattern))$mu
		}else if(gsub("([A-Za-z]+).*", "\\1", drcPattern) == "up"){
			estMu <- ORIClust::up.down(responseData, 
					tapply(X = responseData, INDEX = as.factor(doseData), FUN = mean, na.rm = TRUE), 
					table(doseData), h = parse_number(drcPattern))$mu
		}else if(drcPattern == "complete"){
			estMu <- ORIClust::complete.profile(responseData, 
					tapply(X = responseData, INDEX = as.factor(doseData), FUN = mean, na.rm = TRUE), table(doseData))$mu
		}else{
			stop("The drcPattern should be a valid pattern name.")
		}
		points(as.factor(unique(doseData)), estMu, pch = 13, col = "deeppink", cex = 1.2)
		lines(as.factor(unique(doseData)), estMu, col = "blue4", lwd = 1.5)
	}
}

