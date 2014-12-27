library(shiny)

shinyUI(fixedPage(
	titlePanel("Pancreas Cancer Outcome Predictor"),

	fixedRow(
		column(3,
			wellPanel(
				h4("Curve 1"),
				selectInput("Prognostic1", "Prognostic method:", 
					list("Random Forest" = "rf",
						"Cox" = "cox",
						"Weibull" = "weibull",
						"Kaplan-Meier" = "km")),

				conditionalPanel(
					condition = 'input.Prognostic1 != "km"',
					selectInput("Patient.Sex1", "Patient Sex:",
								list("Male" = "M", "Female" = "F")),

					sliderInput("History.Diagnosis.AgeAt1", "Age at diagnosis:", 
								min = 30, max = 90, value = 60),

					sliderInput("Path.Size1", "Length of longest axis (mm):", 
								min = 5, max = 120, value = 30),

					checkboxInput("Path.LocationHead1", "Tumour in head of pancreas", TRUE)
				),

				checkboxInput("Molec.S100A2.DCThresh1", "S100A2 positivity", FALSE),
				checkboxInput("Molec.S100A4.DCThresh1", "S100A4 positivity", TRUE),

				conditionalPanel(
					condition = 'input.Prognostic1 != "rf"',
					checkboxInput("CIBands1", "Show confidence interval bands", FALSE)
				)
			),
			wellPanel(
				h4("Curve 2"),
				checkboxInput("Show2", "Show curve", FALSE),

				conditionalPanel(
					selectInput("Prognostic2", "Prognostic method:", 
						list("Random Forest" = "rf",
							"Cox" = "cox",
							"Weibull" = "weibull",
							"Kaplan-Meier" = "km")),
					condition = "input.Show2 == true",

					conditionalPanel(
						condition = 'input.Prognostic2 != "km"',
						selectInput("Patient.Sex2", "Patient Sex:",
									list("Male" = "M", "Female" = "F")),

						sliderInput("History.Diagnosis.AgeAt2", "Age at diagnosis:", 
									min = 30, max = 90, value = 60),

						sliderInput("Path.Size2", "Length of longest axis (mm):", 
									min = 5, max = 120, value = 30),

						checkboxInput("Path.LocationHead2", "Tumour in head of pancreas", TRUE)
					),

					checkboxInput("Molec.S100A2.DCThresh2", "S100A2 positivity", FALSE),
					checkboxInput("Molec.S100A4.DCThresh2", "S100A4 positivity", TRUE),

					conditionalPanel(
						condition = 'input.Prognostic2 != "rf"',
						checkboxInput("CIBands2", "Show confidence interval bands", FALSE)
					)
				)
			)
		),
		column(5,
			h4("Survival Curves"),
			fluidRow(
				plotOutput("survPlot")
			)
		),
		column(4,
			h4("Summary Statistics")
		)
	)
))
