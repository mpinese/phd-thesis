library(shiny)

shinyUI(fixedPage(
	titlePanel("Pancreas Cancer Outcome Predictor v1.0"),

	fixedRow(
		column(3,
			wellPanel(
				h4("Curve 1"),
				selectInput("Prognostic1", "Prognostic model:", 
					list("PCOP 1.0" = "pcop10")),

				selectInput("Patient.Sex1", "Patient Sex:",
							list("Male" = "M", "Female" = "F")),

				sliderInput("Path.Size1", "Length of longest axis (mm):", 
							min = 5, max = 120, value = 30),

				checkboxInput("Path.LocationHead1", "Tumour in head of pancreas", TRUE),

				checkboxInput("Molec.S100A2.DCThresh1", "S100A2 positivity", FALSE),
				checkboxInput("Molec.S100A4.DCThresh1", "S100A4 positivity", TRUE),

				checkboxInput("CIBands1", "Show confidence interval bands", FALSE)
			),
			wellPanel(
				h4("Curve 2"),
				checkboxInput("Show2", "Show curve", FALSE),

			 	conditionalPanel(
			 		selectInput("Prognostic2", "Prognostic model:", 
			 					list("PCOP 1.0" = "pcop10")),
			 		condition = "input.Show2 == true",

					selectInput("Patient.Sex2", "Patient Sex:",
								list("Male" = "M", "Female" = "F")),

					sliderInput("Path.Size2", "Length of longest axis (mm):", 
								min = 5, max = 120, value = 30),

					checkboxInput("Path.LocationHead2", "Tumour in head of pancreas", TRUE),

					checkboxInput("Molec.S100A2.DCThresh2", "S100A2 positivity", FALSE),
					checkboxInput("Molec.S100A4.DCThresh2", "S100A4 positivity", TRUE),

					checkboxInput("CIBands2", "Show confidence interval bands", FALSE)
				)
			)
		),
		column(5,
			h4("Survival Curves"),
			fluidRow(
				plotOutput("survPlot")
			)
		)#,
		# column(4,
		# 	h4("Summary Statistics"),
		# 	"TODO"
		# )
	)
))

