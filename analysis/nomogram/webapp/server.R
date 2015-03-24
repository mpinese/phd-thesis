library(shiny)
library(flexsurv)
library(ggplot2)

ENV.PCOP10 = new.env(parent = baseenv())
load("../data/PCOP-1.0.rda", envir = ENV.PCOP10)


library(RColorBrewer)
pal = brewer.pal(3, "Dark2")



calcSurvivalCurve = function(input, which)
{
	if (which == 1)
	{
		newdata = data.frame(
			SexM = input$Patient.Sex1 == "M", 
			LocBody = !input$Path.LocationHead1, 
			SizeCent = input$Path.Size1 - 30, 
			A2 = input$Molec.S100A2.DCThresh1, 
			A4 = input$Molec.S100A4.DCThresh1)
		prognostic = input$Prognostic1
		ci = input$CIBands1
	}
	else
	{
		newdata = data.frame(
			SexM = input$Patient.Sex2 == "M", 
			LocBody = !input$Path.LocationHead2, 
			SizeCent = input$Path.Size2 - 30, 
			A2 = input$Molec.S100A2.DCThresh2, 
			A4 = input$Molec.S100A4.DCThresh2)		
		prognostic = input$Prognostic2
		ci = input$CIBands2
	}

	if (prognostic == "pcop10")
	{
		temp = summary(ENV.PCOP10$fit, newdata = newdata, ci = ci)[[1]]
		if (ci)
		{
			predobj = list(time = temp$time, surv = temp$est, lower = temp$lcl, upper = temp$ucl)
		}
		else
		{
			predobj = list(time = temp$time, surv = temp$est, lower = rep(NA, nrow(temp)), upper = rep(NA, nrow(temp)))
		}
	}

	predobj
}


alphaCol = function(col, alpha)
{
	do.call(rgb, c(as.list(col2rgb(col)/255), alpha = alpha))
}


plotSurvCI = function(times, lower, upper, col)
{
	x = c(times[1], rep(c(times[-1], rev(times)[-1]), each = 2))
	y = c(rep(upper[-length(upper)], each = 2), rep(rev(lower)[-1], each = 2), upper[1])
	polygon(x, y, col = col, border = NA)
}


shinyServer(function(input, output) {
	pred1 <- reactive({
		calcSurvivalCurve(input, 1)
	})

	pred2 <- reactive({
		calcSurvivalCurve(input, 2)
	})

	output$survPlot <- renderPlot({
		plot(pred1()$time/365.25*12, pred1()$surv, type = "s", xlab = "Time from diagnosis (months)", ylab = "Probability of survival", lwd = 2, col = pal[1], xlim = c(0, 60), ylim = c(0, 1), axes = FALSE)
		if (input$CIBands1)
		{
			plotSurvCI(pred1()$time/365.25*12, pred1()$lower, pred1()$upper, alphaCol(pal[1], 0.3))
		}
		box()
		axis(side = 1, at = seq(0, 60, 12))
		axis(side = 2, at = seq(0, 1, 0.1))
		if (input$Show2 == TRUE)
		{
			lines(pred2()$time/365.25*12, pred2()$surv, type = "s", lwd = 2, col = pal[2])
			legend("topright", legend = c("Curve 1", "Curve 2"), col = pal[1:2], lwd = 2, lty = "solid", inset = 0.05)
			if (input$CIBands2)
			{
				plotSurvCI(pred2()$time/365.25*12, pred2()$lower, pred2()$upper, alphaCol(pal[2], 0.3))
			}
		}
		else
		{
			legend("topright", legend = c("Curve 1"), col = pal[1], lwd = 2, lty = "solid", inset = 0.05)
		}
  })
})

