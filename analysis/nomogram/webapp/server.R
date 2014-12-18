library(shiny)
library(randomForestSRC)
library(ggplot2)

load("../03_NSWPCN_subset.rda")
data.all = cbind(time = data.y[,1], event = data.y[,2], data.x.all, data.x.extra)
rsf.molec_preop = rfsrc(
	Surv(time, event) ~ 
		Patient.Sex + 
		History.Diagnosis.AgeAt + 
		Path.LocationBody + 
		Path.Size + 
		Molec.S100A4.DCThresh + 
		Molec.S100A2.DCThresh, 
	data = data.all,
	mtry = 1,
	splitrule = "logrankscore",
	nsplit = 2, 
	ntree = 300)

		# newdata = data.frame(Patient.Sex = factor("M"), History.Diagnosis.AgeAt = 60,
		# 	Path.LocationBody = FALSE, Path.Size = 30, 
		# 	Molec.S100A2.DCThresh = TRUE, Molec.S100A4.DCThresh = FALSE)
		# temp = predict(rsf.molec_preop, newdata)
		# temp$time.interest
		# temp$survival[1,]

shinyServer(function(input, output) {
	pred <- reactive({
		newdata = data.frame(Patient.Sex = factor(input$Patient.Sex, levels = c("M", "F")), History.Diagnosis.AgeAt = input$History.Diagnosis.AgeAt,
			Path.LocationBody = !input$Path.LocationHead, Path.Size = input$Path.Size, 
			Molec.S100A2.DCThresh = input$Molec.S100A2.DCThresh, Molec.S100A4.DCThresh = input$Molec.S100A4.DCThresh)
		predict(rsf.molec_preop, newdata = newdata)
	})

	output$survPlot <- renderPlot({
		plot(pred()$time.interest, pred()$survival, type = "s", xlab = "Time from diagnosis (days)", ylab = "Probability of survival")
  })
})
