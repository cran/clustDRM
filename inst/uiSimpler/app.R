# The ui and server parts of the shiny app
# 
# Author: Vahid Nassiri
###############################################################################


clustDRMapp <- shinyApp(
		ui = 
				navbarPage("clustDRM",
						tabPanel("Data",
								sidebarPanel(
										fileInput("inputData", "Choose input CSV File with the correct format",
												multiple = FALSE,
												accept = c("text/csv",
														"text/comma-separated-values,text/plain",
														".csv")),
										uiOutput("colName"),
										uiOutput("colID"),
										uiOutput("doseID"),
										uiOutput("responseID"),
										uiOutput("addedCovars"),
										uiOutput("plotDataNote"),
										uiOutput("numDigits"),
										uiOutput("isMonotone")),
								mainPanel(
										DT::dataTableOutput("dataTable")
										,
#						textOutput("IDcolp"),
										plotOutput("plotDataOut")
								)),
#						tabPanel("Clustering", sidebarPanel(
#										radioButtons("isMonotone", "Is the dose-response pattern known to be monotone?", 
#												choices  = c("Yes", "No"),
#												selected = "Yes"),
#										selectInput("transform", "Which transform should be applied on the response?", 
#												choices  = c("none", "log","sRoot", "qRoot", "boxcox"),
#												selected = "none"),
#										sliderInput("alpha", "Select the significance level",
#												min = 0.0001, max = 0.5,
#												value = 0.05),
#										actionButton("clustButton", "Perform clustering"),
#										uiOutput("plotClustNote")),
						##				uiOutput("downloadClust")),
#								mainPanel(
#										DT::dataTableOutput("clustRes")
						##,
						##					  textOutput("IDcolp")
#										,
#										plotOutput("plotDataOutClust")
#								)),
						tabPanel("Modelling", sidebarPanel(
										sliderInput("EDp", "Select the effective dose to estimate (EDp)",
												min = 0.1, max = 0.9,
												value = 0.5, step = 0.1),
										selectInput("transform", "Which transform should be applied to the response?", 
												choices  = c("none", "log","sRoot", "qRoot", "boxcox"),
												selected = "none"),
										uiOutput("selThreshold"),
										actionButton("fitButton", "Estimate EDp"),
										uiOutput("downloadRes"),
										uiOutput("plotFitNote")	
								),
								mainPanel(
										DT::dataTableOutput("fitRes")
#,
#					  textOutput("IDcolp")
										,
										plotOutput("plotDatafit")
								)),
						tabPanel("Instructions",
								tags$iframe(style = "width: 100%; border: none; height: 900px;", src = "helpApp.html"))
#						tabPanel("Simulation",
#								sidebarPanel(
#										uiOutput("pilotID"),
#										textInput("simDoseLevels", "Enter dose levels (comma delimited)", "0,5,10"),
#										textInput("simNumRep", "Enter number of replications per dose level (comma delimited)", "24,3,3"),
#										numericInput('numSim', "Enter number of replications", value = 100, min = 1),
#										numericInput('standardDeviation', "Enter standard deviation of the response", 
#												value = 0.1, min = 0.00001, step = 0.05),
#										sliderInput("EDpSim", "Select the effective dose to be estimated (EDp)",
#												min = 0.1, max = 0.9,
#												value = 0.5, step = 0.1),
#										checkboxGroupInput("funcList", "Select the models to be used in the simulation study",
#												c("linear", "linlog", "exponential", "emax", "sigEmax", "logistic", "betaMod","quadratic"), 
#												selected = c("linear", "linlog", "exponential", "emax", "sigEmax", "logistic", "betaMod","quadratic")),
#										actionButton("simButton", "Run simulation"),
#										uiOutput("simPlot")
#								),
#								mainPanel(
#										plotOutput("plotSimRes", height = "600px")
				##						,
				##										textOutput("printTest")
#								))
				
				## For modelling we need two things: p and the threshold. Threshold is optional.
				## For plot we shoud let the user select from all comoounds. perhaps we give it a list
				## as a vector with the compound names to select, etc.
				)
		
		,
		server = function(input, output) {
			round_df <-
			function(df, digits = 3) {
				nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
				
				df[,nums] <- round(df[,nums], digits = digits)
				
				(df)
			}
			inputDataActive <- reactive({
						read.csv(input$inputData$datapath)
					})	
			## Selecting name, ID, dose and response columns
			
			output$colName <- renderUI({
#						if(is.null(input$inputData))
#							return()
#						
						req(input$inputData)
						selectInput("colName", "Select the name", 
								choices  = colnames(inputDataActive()),
								selected = NULL)
					})
			
			nameCol <- reactive({
						which(colnames(inputDataActive()) == input$colName)
					})	
			
			output$colID <- renderUI({
#						if(is.null(input$inputData))
#							return()
						req(input$inputData)
						
						selectInput("colID", "Select the ID", 
								choices  = colnames(inputDataActive())[-nameCol()],
								selected = NULL)
					})
			
			IDcol <- reactive({
						which(colnames(inputDataActive()) == input$colID)
					})	
			
			output$doseID <- renderUI({
#						if(is.null(input$inputData))
#							return()
#						
						req(input$inputData)
						
						selectInput("doseID", "Select the dose", 
								choices  = colnames(inputDataActive())[-c(nameCol(),IDcol())],
								selected = NULL)
					})
			
			doseCol <- reactive({
						which(colnames(inputDataActive()) == input$doseID)
					})	
			
			output$responseID <- renderUI({
#						if(is.null(input$inputData))
#							return()
						req(input$inputData)
						
						
						selectInput("responseID", "Select the response", 
								choices  = colnames(inputDataActive())[-c(nameCol(),IDcol(), doseCol())],
								selected = NULL)
					})
			
			responseCol <- reactive({
						which(colnames(inputDataActive()) == input$responseID)
					})	
			
			
			output$addedCovars <- renderUI({
#						if(is.null(input$inputData))
#							return()
						req(input$inputData)
						
						selectInput("addedCovars", "Select the added covariates", 
								choices  = colnames(inputDataActive())[-c(nameCol(),IDcol(), doseCol(), responseCol())],
								selected = NULL, multiple = TRUE)
					})
			
			
			output$numDigits <-renderUI({
#						if(is.null(input$inputData))
#							return()
#						
						req(input$inputData)
						
						sliderInput("numDigits", "Select the number of digits to be displayed", 0, 15, 2)
					})
			
			output$isMonotone <- renderUI({
#						if(is.null(input$inputData))
#							return()
#						
						req(input$inputData)
						
						radioButtons("isMonotone", "Is the dose-response pattern known to be monotone?", 
								choices  = c("Yes", "No"),
								selected = "Yes")
					})
			
			
#			output$transform <- renderUI({
#						if(is.null(input$inputData))
#							return()
#						
#						selectInput("transform", "Which transform should be applied on the response?", 
#								choices  = c("none", "log","sRoot", "qRoot", "boxcox"),
#								selected = "none")
#					})
			
			
			
			output$dataTable <- DT::renderDataTable({
						req(input$inputData)
						DT::datatable(round_df(inputDataActive(), input$numDigits), selection = "single")
					})
			
			
			dataIDToPlot <- reactive({
						req(input$dataTable_rows_selected)
						inputDataActive()[input$dataTable_rows_selected, IDcol()]
					})
			
			
			output$plotDataNote <- renderUI({
						# if(is.null(input$colID) | is.null(input$doseID) |is.null(input$responseID) )
						# 	return()
						
						req(input$colID)
						req(input$doseID)
						req(input$responseID)
						
						p("Select a subject to plot the dose-response relation with average response (shown using a cross) per dose (tiny differences between plots is due to generated jitter).")
					})
			
			addCovars <- reactive({
						if (is.null(input$addedCovars)){
							as.formula("~1")
						}else{
							as.formula(paste("~", paste(input$addedCovars, collapse = "+")))
						}
					})
			
			output$plotDataOut <- renderPlot({
						req(input$inputData, input$colID, input$doseID, input$responseID, input$dataTable_rows_selected)
						plotDoseResponseData(inputDataActive(), doseCol(), responseCol(), IDcol(), dataIDToPlot())
					})
			
			isMonotone <- reactive({input$isMonotone})
			transform <- reactive({input$transform})
			toInput <- reactive({inputDataMaker(doseCol(), responseCol(), IDcol(), inputDataActive())})
			
#			
#			toDisplay <- eventReactive(input$clustButton, {
#		
#					})
			
			
#	if (!is.null(clusteringResults())){
#		colnames(clusteringResults()) <- c("ID", "Identified pattern", "Adjusted p-value")
#	}
#	
			
#			output$clustRes <- DT::renderDataTable({
#						req(toDisplay())
#						DT::datatable(round_df(toDisplay(), input$numDigits), selection = "single")
#					})	
#			
#			dataIDToPlotClust <- reactive({
#						req(input$clustRes_rows_selected)
#						toDisplay()[input$clustRes_rows_selected, 2]
#					})
#			
#			dataPatternToPlotClust <- reactive({
#						req(input$clustRes_rows_selected)
#						toDisplay()[input$clustRes_rows_selected, 3]
#					})
			
			
#			output$plotClustNote <- renderUI({
#						if(is.null(toDisplay()))
#							return()
#						p("Select a subject to plot the dose-response relation with estimated average per dose based on the identified pattern (tiny differences between plots is due to generated jitter).")
#					})
#			
#			output$plotDataOutClust <- renderPlot({
#						req(input$inputData, input$colID, input$doseID, input$responseID, input$clustRes_rows_selected)
#						plotDoseResponseData(inputDataActive(), doseCol(), responseCol(), IDcol(), subjectID = dataIDToPlotClust(), 
#								addMean = FALSE, drcPattern = dataPatternToPlotClust())
#					})
#			
#			
			
			
			output$selThreshold <- renderUI({
#						if(is.null(input$inputData))
#							return()
						req(input$inputData)
						
						numericInput("selThreshold", "Enter the toxic dose threshold", 
								value = NA,
								min = min(inputDataActive()[, doseCol()[]]), max = max(inputDataActive()[, doseCol()[]]))
					})
			
			
			
			
			
			## making what should be displayed after modelling
			
			toDisplayMod <- eventReactive(input$fitButton, {
						## clustering
						withProgress(message = "Clustering is in progress",
								detail = "depending on the number of subjects, it may take several minutes...", value = 0, {
									if (isMonotone() == "Yes"){
										clusteringResults0 <- monotonePatternClustering (inputData = toInput()$inputData, colsData = toInput()$colsData ,
												colID = toInput()$colID, 
												doseLevels = toInput()$doseLevels, numReplications = toInput()$numReplicates, transform = transform(), 
												BHorBY = TRUE, SAM = FALSE, testType = "E2",
												adjustType = "BH", FDRvalue = c(0.05, 0.05), nPermute= c(1000, 1000), fudgeSAM = "pooled",
												useSeed = c(NULL, NULL), theLeastNumberOfTests = 1, na.rm = FALSE, imputationMethod = "mean")
										clusteringResults1 <- data.frame(clusteringResults0$resultsBH$E2, 
												clusteringResults0$subjectsPatterns[clusteringResults0$subjectsPatterns!= "flat"])
										
										idSel <- clusteringResults0$selectedSubjectsBH$CompID[clusteringResults0$subjectsPatterns!= "flat"]
										
										if (length(idSel) > 0){
											namesComp <- factor(rep(NA, length(idSel)), levels = levels(inputDataActive()[,nameCol()]))
											for (iID in seq_along(idSel)){
												namesComp[iID] <- inputDataActive()[inputDataActive()[,IDcol()] == idSel[iID], nameCol()][1]
											}
										}else{
											stop("None of the compounds have a non-flat pattern!")
										}
										
										clusteringResults2 <- data.frame(namesComp,
												clusteringResults1[,c(2, 5, 4)])
										colnames(clusteringResults2) <- c("Name", "ID", "Identified pattern", "Adjusted p-value")
										clusteringResults2
									}else{
										clusteringResults0 <- generalPatternClustering(inputData = toInput()$inputData, colsData = toInput()$colsData ,colID = toInput()$colID , 
												doseLevels = toInput()$doseLevels, numReplications = toInput()$numReplicates, na.rm = FALSE, imputationMethod = "mean",
												ORICC = "two", transform = "none",plotFormat = "eps", LRT = FALSE, MCT = TRUE,
												adjustMethod = "BH",
												nPermute = 1000, useSeed = NULL, theLeastNumberOfMethods = 2, alpha = 0.05, 
												nCores = parallel::detectCores(all.tests = FALSE, logical = TRUE)-1)
										clusteringResults1 <- clusteringResults0$clusteringORICC2Results$resultsMCT[clusteringResults0$clusteringORICC2Results$resultsMCT[,1]!= "flat",
												-c(3,6)]
										idSel <- clusteringResults0$clusteringORICC2Results$resultsMCT$ID[clusteringResults0$clusteringORICC2Results$resultsMCT[,1]!= "flat"]
										if (length(idSel) > 0){
											namesComp <- factor(rep(NA, length(idSel)), levels = levels(inputDataActive()[,nameCol()]))
											for (iID in seq_along(idSel)){
												namesComp[iID] <- inputDataActive()[inputDataActive()[,IDcol()] == idSel[iID], nameCol()][1]
											}
										}else{
											stop("None of the compounds have a non-flat pattern!")
										}
										clusteringResults2 <- data.frame(namesComp,
												clusteringResults1[,c(2, 1, 4)])
										colnames(clusteringResults2) <- c("Name", "ID", "Identified pattern", "Adjusted p-value")
										
									}
									
								})
						## Modelling
						withProgress(message = "Fitting model is in progress",
								detail = "depending on the number of subjects, it may take several minutes...", value = 0, {
									fittedMod <-fitDRM (inputDataActive(), doseCol(), responseCol(), IDcol(), subsettingID = clusteringResults2[,2], 
											transform = transform(), addCovars = addCovars(), patternClusters = clusteringResults2[,3], 
											EDp = input$EDp, addCovarsVar = TRUE, alpha = 0.05, na.rm = FALSE, imputationMethod = "mean", 
											nCores = parallel::detectCores(all.tests = FALSE, logical = TRUE)-1)
									
									estimatedEDp <- cbind(as.numeric(c(as.character(fittedMod$estEDpMonotone$CompID), as.character(fittedMod$estEDpNonmonotone$CompID))),
											c(fittedMod$estEDpMonotone$modelAveragingAIC, fittedMod$estEDpNonmonotone$modelAveragingAIC))
									
									idSel <- estimatedEDp[,1]
									if (length(idSel) > 0){
										namesComp <- factor(rep(NA, length(idSel)), levels = levels(inputDataActive()[,nameCol()]))
										for (iID in seq_along(idSel)){
											namesComp[iID] <- inputDataActive()[inputDataActive()[,IDcol()] == idSel[iID], nameCol()][1]
										}
									}else{
										stop("Model fitting is not possible! None of the compounds have a monotone or umbrella-shape pattern.")
									}
									
									if (is.na(input$selThreshold)){
										estEDp0 <- data.frame(namesComp,
												estimatedEDp[order(estimatedEDp[,1]),])
										
										
										colnames(estEDp0) <- c("Name", "ID", "estimated EDp")
										estEDp0
									}else{
										isToxic0 <- estimatedEDp[order(estimatedEDp[,1]),2] < input$selThreshold
										isToxic <- rep("No", length(isToxic0))
										isToxic[isToxic0] <- "Yes"
										estEDp0 <- data.frame(namesComp,
												estimatedEDp[order(estimatedEDp[,1]),], 
												isToxic)
										
										
										colnames(estEDp0) <- c("Name", "ID", "estimated EDp", "is toxic?")
										estEDp0
									}
									
									
								})
					})
			
			
			output$fitRes <- DT::renderDataTable({
						req(toDisplayMod())
						DT::datatable(round_df(toDisplayMod(), input$numDigits), selection = "single")
					})	
			
			output$downloadRes <- renderUI({
						req(toDisplayMod())
						downloadButton('estResDown', 'Download results')
						
					})
			
			output$estResDown <- downloadHandler(
					filename = function() {
						paste(input$responseID,"_results_", Sys.time(),".csv", sep = "")
					},
					content = function(file) {
						write.csv(toDisplayMod(), file, row.names = FALSE)
					}
			)
			
			
			dataIDToPlotFit <- reactive({
						req(input$fitRes_rows_selected)
						toDisplayMod()[input$fitRes_rows_selected, 2]
					})
			
			EDpToPlot <- reactive({
						req(input$fitRes_rows_selected)
						toDisplayMod()[input$fitRes_rows_selected, 3]
					})
			
			## given the doses, here we define in which interval the estimated
			## EDp would fall
			
			uniqueDoses <- reactive(unique(inputDataActive()[,doseCol()]))
			
			doseIntervals <- reactive(c(which(uniqueDoses() > EDpToPlot())[1]-1, 
							which(uniqueDoses() > EDpToPlot())[1]))
			
			
			output$plotFitNote <- renderUI({
#						if(is.null(toDisplayMod()))
#							return()
#						
						req(toDisplayMod())
						
						p("Select a subject to plot the dose-response relation with estimated EDp -red dashed line- (tiny differences between plots is due to generated jitter).")
					})
			
			output$plotDatafit <- renderPlot({
						req(input$inputData, input$colID, input$doseID, input$responseID, input$fitRes_rows_selected)
						plotDoseResponseData(inputDataActive(), doseCol(), responseCol(), IDcol(), subjectID = dataIDToPlotFit(), 
								addMean = FALSE, drcPattern = NULL)
						abline(v = doseIntervals()[1] + ((EDpToPlot() - uniqueDoses()[doseIntervals()[1]])/ 
											(uniqueDoses()[doseIntervals()[2]]-uniqueDoses()[doseIntervals()[1]])), lty = 2, col = "red")
					})
			
			
			
			
			
#			## SIMULATINS
#			
#			output$pilotID <- renderUI({
#						if(is.null(input$inputData))
#							return()
#						
#						selectInput("pilotID", "Select a subject to be used as pilot dataset", 
#								choices  = unique(inputDataActive()[, nameCol()]),
#								selected = NULL, multiple = FALSE)
#					})
#			
#			pilotData <- reactive({
#						req(input$pilotID)
#						inputDataActive()[which(inputDataActive()[, nameCol()] == input$pilotID),c(doseCol(), responseCol())]
#					})
#			
#			doseLevels <- reactive({as.numeric(unlist(strsplit(input$simDoseLevels,",")))})
#			numReplications <- reactive({as.numeric(unlist(strsplit(input$simNumRep,",")))})
#			
#			
#			
#			
#			
#			simRes <- eventReactive(input$simButton, {
#						withProgress(message = "Simulation is in progress",
#								detail = "it may take several minutes...", value = 0, {
#									
#									simulEvalDRM(pilotData(), doseLevels(), numReplications(), numSim = input$numSim, 
#											standardDeviation = input$standardDeviation, EDp = input$EDpSim,
#											funcList = input$funcList)
#									
#								})
#					})
#			
#			output$simPlot <- renderUI({
#						if(is.null(simRes()))
#							return()			
#						radioButtons("quantity2Plot", "Select a measure to summrize simulation results",
#								choices = c("bias" = "bias",
#										"MSE" = "mse",
#										"relatrive Bias" = "relativeBias",
#										"absolute bias" = "absBias",
#										"absolute relative bias" = "absRelativeBias"),
#								selected = "mse")
#					})
#			
#			
#			output$plotSimRes <- renderPlot({
#						req(input$quantity2Plot)
#						simToPlot <- simRes()
#						plotSimulDRM(simToPlot, quantity2Plot = input$quantity2Plot)
#					})
#			
			
#			output$printTest <- renderPrint({
#						req(input$quantity2Plot)
#						simToPlot <- simRes()
#						input$quantity2Plot
			
#					})
			
			
		}
)
