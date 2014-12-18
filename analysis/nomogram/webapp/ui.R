library(shiny)

shinyUI(
pageWithSidebar(
	headerPanel("Survival Test"),

	sidebarPanel(
		selectInput("Patient.Sex", "Patient Sex:",
					list("Male" = "M", "Female" = "F")),

		sliderInput("History.Diagnosis.AgeAt", "Age at diagnosis:", 
					min = 30, max = 90, value = 60),

		sliderInput("Path.Size", "Length of longest axis (mm):", 
					min = 5, max = 120, value = 30),

		checkboxInput("Path.LocationHead", "Tumour in head of pancreas", TRUE),

		checkboxInput("Molec.S100A2.DCThresh", "S100A2 positivity", FALSE),
		checkboxInput("Molec.S100A4.DCThresh", "S100A4 positivity", TRUE)
    ),

	mainPanel(
    	plotOutput("survPlot")
	)
))
