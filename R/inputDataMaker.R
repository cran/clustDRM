# TODO: Add comment
# 
# Author: Vahid Nassiri
###############################################################################


#' Creating suitable inputData for clustering of the dose-response curve patterns
#' 
#' function to create needed information as the input of the functions to cluster 
#' dose-response cruve patterns.
#' 
#' @param dose either a single string or a scalar, indicating the name of the dose column or its index.
#' @param response either a single string or a scalar, indicating the name of the response column or its index.
#' @param ID either a single string or a scalar, indicating the name of the ID column or its index.
#' @param inputDataset a data frame containing the input dataset, it should at least include dose, response, and ID
#' 
#' @return a list with the following elements:
#' 
#' inputDataset: includes the ID (first column), and the response for all doses with their replications for each subject as rows.
#' doseLevels: unique dose levels
#' numReplicatrions: number of replicatios per each unique dose level.
#' colsData: the index of columns with responses.
#' colID: the index of ID column.
#' 
#' @details Note that the output of this function can be feed into the function for clustering dose-response curve patterns.
#'
#' @importFrom MCPMod genDFdata 
#' @examples  
#' ## gnerating data
#' set.seed(11)
#' doses2Use <-  c(0, 5, 20)
#' numRep2Use <- c(3, 3, 3)
#' generatedData <- cbind(rep(1,sum(numRep2Use)), 
#' MCPMod::genDFdata("logistic",c(5, 3, 10, 0.05), 
#' doses2Use, numRep2Use, 1), 
#' 		matrix(rnorm(1*sum(numRep2Use)), sum(numRep2Use), 1))
#' colnames(generatedData) <- c("ID", "dose", "response", "x1")
#' for (iGen in 2:15){
#' 	genData0 <- cbind(rep(iGen,sum(numRep2Use)), 
#' MCPMod::genDFdata("logistic",c(5, 3, 10, 0.05), 
#' doses2Use, numRep2Use, 1), 
#' 			matrix(rnorm(1*sum(numRep2Use)), sum(numRep2Use), 1))
#' 	colnames(genData0) <- c("ID", "dose", "response", "x1")
#' 	generatedData <- rbind(generatedData, genData0)
#' }
#' ## transforming it for clustering
#' toInput <- inputDataMaker(2, 3, 1, generatedData)
#' 
#' @author Vahid Nassiri, and Yimer Wasihun
#' @export
inputDataMaker <- function(dose, response, ID, inputDataset){
	## split the input dataframe into a list of daraframes each with the data regarding each subject
	if (is.character(ID)){
		colID <- which(names(inputDataset) == ID)
	}else{
		colID <- ID
	}
	if (is.character(dose)){
		doseCol <- which(names(inputDataset) == dose)
	}else{
		doseCol <- dose
	}
	if (is.character(response)){
		responseCol <- which(names(inputDataset) == response)
	}else{
		responseCol <- response
	}
	inputDataList <- split(inputDataset, f = unlist(inputDataset[,colID]))
	## find the doseLevels and numReplications. Note that here we assume 
	## different compounds have the same number of replications, as well as 
	doseLevels <- unique(inputDataList[[1]][,doseCol])
	numReplicates <- table(inputDataList[[1]][,doseCol])
	inputData0 <- matrix(NA, length(inputDataList), nrow(inputDataList[[1]]))
	for (iSubj in seq_along(inputDataList)){
		inputData0[iSubj, ] <- inputDataList[[iSubj]][, responseCol]
	}
	inputData <- data.frame(unique(inputDataset[,colID]), inputData0)
	colsData <- 2:ncol(inputData)
	doseRep <- cbind(inputDataList[[1]][,doseCol], unlist(lapply(numReplicates, seq)))
	names(inputData) <- c(names(inputDataset)[colID], paste(doseRep[,1], doseRep[,2], sep = ":"))
	return(list(inputData = inputData, doseLevels = doseLevels, numReplicates = numReplicates, colsData = colsData, colID = 1))
}

