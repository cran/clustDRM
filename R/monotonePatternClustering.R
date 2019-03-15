# TODO: Add comment
# 
# Author: Vahid Nassiri
###############################################################################

#' clustering dose-response curves based on their pattern when it is 
#' known to be monotone.
#' function to cluster dose-response curves based on their pattern.
#' 
#' @param inputData data matrix which should incluide ID's of the subjects, as well as the measurements (gene expressions, etc.) 
#' for all replications of different as columns. 
#' @param colsData vector indicating the idex of columns in the inputData which correspond to the measurement for different replications 
#' of different doses.
#' @param colID scalar indicating the index of column corresponding to data ID.
#' @param doseLevels vector with dose levels.
#' @param transform single string indicating what kind of transform should be applied on the response data. 
#' It takes "none" (no transform, dafault), "log" (natural log), "sRoot (square root), and "qRoot" (cubic root), and 
#' "boxcox" (Box-Cox transformation).
#' @param numReplications vector wit hthe same length as doseLevels with number of replications for each dose. 
#' @param BHorBY logical inidicating whether monotonicity tests (specified in argument testType) using BH or BY modifications 
#' should be performed. Default is TRUE.
#' @param SAM logical indicating whether a SAM procedure should be perfmored. Default is FALSE.
#' @param testType string a subset of c("E2", "Williams", "Marcus", "M", "ModifM"), indicating the monotonicity tests 
#' which should be applied.
#' @param adjustType method of adjustment for multi-plicity in case BHorBY = TRUE. It takes values "BH" abd "BY" with "BH" as the
#' default.
#' @param FDRvalue a numerical vector of length 2 indicating the FDR values for BHorBY and SAM, the default is 0.05 for both. 
#' @param nPermute a numerical vector of length 2 indicating number of permutatation for BHorBY and SAM, the default for 
#' both is 1000.
#' @param fudgeSAM single string takes value from ("pooled", "none") specified the fudge factor in SAM test statistic. 
#' The default is "pooled".
#' @param useSeed a vector of lkength two specifying the seed value for BHorBY and SAM, the default is NULL for both.
#' @param theLeastNumberOfTests A scalar indicating the minimum number of tests which should approve a monotone trend to 
#' consider a trend montone. The default is 5, i.e., all of the tests should agree on the monotonicity.
#' @param na.rm logical variable indicatign whether missing values should be removed (TRUE) or not (FALSE, default)
#' @param imputationMethod signle string taking calues from "mean" (default), and "median", which indicates how the missing values should be
#' treated. "mean" would replace them with the mean of the observed ones, and "median" will use median of them for imputation.
#' 
#' @return a list with the following objects:
#' 
#' selectedSubjects: provides the ID and indentified patterns for the subjects which are selected based on the results of 
#' various tests and theLeastNumberOfTests.
#' 
#' subjectsPatterns: a vector of the same length as the number of subjects in the input dataset which indicates the identified
#' patterns for all subjects (including flat ones).
#' 
#' resultsBH: a list with the results of selected tests (if BHorBY = TRUE, NULL otherwise).
#' 
#' resultsSAM: a list with results of SAM procedure (if SAM = TRUE, NULL otherwise).
#' 
#' selectedSubjectsBH: a data frame of all of the subjects with then number of tests select them based on adjusted BH or BY methods.
#' 
#' selectedSubjectsSAM: a data frame of all of the subjects with then number of tests select them based on SAM procedure
#'  
#' @examples
#' ## gnerating data, a sample of size 20
#' set.seed(11)
#' doses2Use <-  c(0, 5, 20)
#' numRep2Use <- c(3, 3, 3)
#' generatedData <- cbind(rep(1,sum(numRep2Use)), 
#' MCPMod::genDFdata("logistic",c(5, 3, 10, 0.05), 
#' doses2Use, numRep2Use, 1), 
#' 		matrix(rnorm(1*sum(numRep2Use)), sum(numRep2Use), 1))
#' colnames(generatedData) <- c("ID", "dose", "response", "x1")
#' for (iGen in 2:20){
#' 	genData0 <- cbind(rep(iGen,sum(numRep2Use)), 
#' MCPMod::genDFdata("logistic",c(5, 3, 10, 0.05), 
#' doses2Use, numRep2Use, 1), 
#' 			matrix(rnorm(1*sum(numRep2Use)), sum(numRep2Use), 1))
#' 	colnames(genData0) <- c("ID", "dose", "response", "x1")
#' 	generatedData <- rbind(generatedData, genData0)
#' }
#' ## transforming it for clustering
#' toInput <- inputDataMaker(2, 3, 1, generatedData)
#' ## monotone pattern clustering
#' monotonePatternClust <- monotonePatternClustering (inputData = 
#' toInput$inputData, colsData = toInput$colsData ,
#' 		colID = toInput$colID, doseLevels = toInput$doseLevels, 
#' numReplications = toInput$numReplicates, 
#' 		BHorBY = TRUE, SAM = FALSE, testType = c("E2"),
#' 		adjustType = "BH", FDRvalue = c(0.05, 0.05), 
#' nPermute= c(100, 100), fudgeSAM = "pooled",
#' 		useSeed = c(NULL, NULL), theLeastNumberOfTests = 1, 
#' na.rm = FALSE, imputationMethod = "mean")
#' 
#' @seealso \href{https://www.rdocumentation.org/packages/IsoGene/versions/1.0-24/topics/IsoTestBH}{IsoGene}
#' \href{https://www.rdocumentation.org/packages/IsoGene/versions/1.0-24/topics/IsoTestSAM}{IsoGene}
#' \href{https://www.rdocumentation.org/packages/ORCME/versions/2.0.2/topics/monotoneDirection}{ORCME}
#' 
#' @import IsoGene
#' @import ORCME
#' @importFrom MCPMod genDFdata 
#' @author Vahid Nassiri, and Yimer Wasihun.
#' @export
monotonePatternClustering <- function(inputData, colsData ,colID, 
		doseLevels, numReplications, transform = c("none", "log","sRoot", "qRoot", "boxcox"), BHorBY = TRUE, SAM = FALSE, testType = c("E2", "Williams", "Marcus", "M", "ModifM"),
		adjustType = c("BH", "BY"), FDRvalue = c(0.05, 0.05), nPermute= c(1000, 1000), fudgeSAM = c("pooled", "none"),
		useSeed = c(NULL, NULL), theLeastNumberOfTests = 5, na.rm = FALSE, imputationMethod = c("mean", "median")){
	
	adjustType <- match.arg(adjustType)
	## keeping original input data
	origInputData <- inputData
	## As we are constructing that outseveles, so we know it's in column 3
	## finding subject with missing values
	idxNA <- which(is.na(inputData), arr.ind=TRUE)
	idxNARows <- unique(idxNA[,1])
	## either remove subjects with missing values or replace them with
	## the mean or median of the avilable values.
	if (na.rm){
		inputData <- inputData[-idxNARows,]
	}else{
		missingSubjects <- inputData[idxNARows,-colID]
		imputedValues <- apply(missingSubjects, 1, get(imputationMethod), na.rm = TRUE)
		for (iNA in 1:length(idxNARows)){
			inputData[idxNARows[iNA], idxNA[idxNA[,1]== idxNARows[iNA],2]] <- imputedValues[iNA]
		}
	}		
	## for the ones which are vectors, the first argument is for BH, second for SAM, so if any of them is set to TRUE,
	## it's corresponding arguments should also be specified.
	
	## theLeastNumberOfTests cannot be larger than length(testType) * sum(BHorBY + SAM) and it indicates the 
	## least number of test which should select a compound to let it in. If it is set to 10, that means all of the methods
	## should have selected that subject, it's equal to what is in the code now.
	testType <- match.arg(testType, c("E2", "Williams", "Marcus", "M", "ModifM"), several.ok = TRUE)
	if (theLeastNumberOfTests <1 | theLeastNumberOfTests > (length(testType) * sum(BHorBY + SAM))){
		stop(paste("theLeastNumberOfTests should be at least 1 and at most ", 
						length(testType) * sum(BHorBY + SAM)))
	}
	transform <- match.arg(transform, c("none", "log","sRoot", "qRoot", "boxcox"))
	## specifying parameters
	FDRvalueBH <- FDRvalue[1]
	FDRvalueSAM <- FDRvalue[2]
	nPermuteBH <- nPermute[1]
	nPermuteSAM <- nPermute[2]
	useSeedBH <- useSeed[1]
	useSeedSAM <- useSeed[2]
	## transforming data if needed
	if (transform == "none"){
		inputData[,colsData] <- inputData[,colsData]
	}else if(transform == "log"){
		if (any(unlist(inputData[,colsData]) == 0 )){
			inputData[,colsData] <- log(inputData[,colsData] + 1)
			warning("Due to zero responses, log(x+1) is used instead of log.")
		}
		if(any(unlist(inputData[,colsData]) <0)){
			stop("log transofrm is not possible with negative responses.")
		}
		if(all(unlist(inputData[,colsData])>0)){
			inputData[,colsData] <- log(inputData[,colsData])
		}
	}else if (transform == "sRoot"){
		if (any(unlist(inputData[,colsData])<0)){
			stop("square root is not possible with negative responses.")
		}else{
			inputData[,colsData] <- sqrt(inputData[,colsData])
		}
	}else if (transform == "qRoot"){
		inputData[,colsData] <- sign(inputData[,colsData]) * (abs(inputData[,colsData])^(1/3))
	}else if(transform  == "boxcox"){
		inputData[,colsData] <- apply(inputData[,colsData], 2, boxCoxTransform)
	}
	
	
	
	## repeating dose levels for their number of replications
	dose <- rep(doseLevels, numReplications)
	subjectsID <- inputData[,colID]
	finalSelectedSubjects <- subjectsID
	## pre-defining outputs
	resultsBH <- NULL
	selectedSubjectsBH <- NULL
	resultsSAM <- NULL
	selectedSubjectsSAM <- NULL
	if (BHorBY){
		if (!is.null(useSeedBH)){
			set.seed(useSeedBH)
		}
		rawPvalue <- IsoGene::IsoRawp(dose, as.data.frame(inputData[,colsData]), niter = nPermuteBH, progressBar = FALSE)
		## we are only interested in te two sided p-values
		twoSidedPvalue <- rawPvalue[[2]]
		names(twoSidedPvalue)[1] <- names(inputData)[colID]
		## Extract the compounds significantly monotonic using all the tests
		resultsBH <- mapply(IsoGene::IsoTestBH, rp = rep(list(twoSidedPvalue), length(testType)), 
				FDR = rep(list(FDRvalueBH), length(testType)), 
				type = rep(list(adjustType), length(testType)),
				stat = as.list(testType), SIMPLIFY = FALSE)
		names(resultsBH) <- testType
		selectedSubjects <- matrix(0,nrow(inputData), length(testType))
		for (iSelect in seq_along(testType)){
			selectedSubjects[,iSelect] <- subjectsID %in% resultsBH[[iSelect]]$Probe.ID
		}
		selectedSubjectsBH <- data.frame(subjectsID, selectedSubjects, apply(selectedSubjects, 1, sum))
		names(selectedSubjectsBH) <- c(colnames(inputData)[colID], testType, "selectedByHowmanyTests")
		finalSelectedSubjects <- cbind(finalSelectedSubjects, selectedSubjects)
	}
	if (SAM){
		## SAM
		if (!is.null(useSeedSAM)){
			set.seed(useSeedSAM)
		}
		resultsSAM <- mapply(IsoGene::IsoTestSAM, 
				x = rep(list(dose), length(testType)), 
				y = rep(list(inputData [, colsData]), length(testType)), 
				fudge = rep(list(fudgeSAM), length(testType)), 
				niter = rep(list(nPermuteSAM), length(testType)), 
				FDR = rep(list(FDRvalueSAM), length(testType)), 
				stat = as.list(testType), SIMPLIFY = FALSE)
		
		names(resultsSAM) <- resultsSAM
		selectedSubjects <- matrix(0,nrow(inputData), length(testType))
		for (iSelect in seq_along(testType)){
			selectedSubjects[,iSelect] <- subjectsID %in% resultsSAM[[iSelect]][[1]]$Probe.ID
		}
		selectedSubjectsSAM <- data.frame(subjectsID, selectedSubjects, apply(selectedSubjects, 1, sum))
		names(selectedSubjectsSAM) <- c(colnames(inputData)[colID], testType, "selectedByHowmanyTests")
		finalSelectedSubjects <- cbind(finalSelectedSubjects, selectedSubjects)
	}
	## Selecting the subject sill then go to testing for increasing or decreasing trends
	if (!SAM & !BHorBY){
		subjectsToUse <- subjectsID
	}else{
		subjectsToUse <- subjectsID[which(apply(as.matrix(finalSelectedSubjects[,-1]),1, sum) >= theLeastNumberOfTests)]
	}
	## testing to see if the selected subjects are increasing or decreasing
	if (length(theLeastNumberOfTests) > 0){
		directionOfPattern <- ORCME::monotoneDirection(geneData = inputData[subjectsToUse, colsData],
				doseData = dose)
		FinalDirectionOfPattern <- directionOfPattern$direction
		FinalDirectionOfPattern[FinalDirectionOfPattern=="dn"] <- "decreasing"
		FinalDirectionOfPattern[FinalDirectionOfPattern=="up"] <- "increasing"
		## the final selected subjects with their patterns
		selectedSubjects <- data.frame(subjectsToUse, FinalDirectionOfPattern)
		colnames(selectedSubjects) <- c(names(inputData)[colID], "pattern")
		## determining the direction, also change it from dn and up to
		## decreasing and increasing
		subjectsPatterns <- rep("flat", nrow(inputData))
		subjectsPatterns[subjectsToUse] <- FinalDirectionOfPattern
	}else{
		warning("No subject is selected to be used!")
	}
	toReturn <- list(selectedSubjects = selectedSubjects, subjectsPatterns = subjectsPatterns, resultsBH = resultsBH, 
			selectedSubjectsBH = selectedSubjectsBH, resultsSAM = resultsSAM,
			selectedSubjectsSAM = selectedSubjectsSAM)
	return(toReturn)
}

