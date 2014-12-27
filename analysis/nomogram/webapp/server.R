library(shiny)
library(randomForestSRC)
library(survival)
library(ggplot2)
library(rms)

load("../03_NSWPCN_subset.rda")
data.all = cbind(time = data.y[,1], event = data.y[,2], data.x.all, data.x.extra)

dd <<- datadist(data.all)
options(datadist = "dd")

rsf.molec_preop = rfsrc(
	Surv(time, event) ~ 
		Patient.Sex + 
		History.Diagnosis.AgeAt.Cent + 
		Path.LocationBody + 
		Path.Size.Cent + 
		Molec.S100A4.DCThresh + 
		Molec.S100A2.DCThresh, 
	data = data.all,
	mtry = 1,
	splitrule = "logrankscore",
	nsplit = 2, 
	ntree = 1000)

km.molec_preop = survfit(
	Surv(time, event) ~ 
		Molec.S100A4.DCThresh + 
		Molec.S100A2.DCThresh, 
	data = data.all)
km.molec_preop.strata = rep(gsub(" ", "", names(km.molec_preop$strata)), km.molec_preop$strata)

# cox.molec_preop.null = cph(
# 	Surv(time, event) ~ 1,
# 	data = data.all, x = TRUE, y = TRUE)
# par(ask = TRUE)
# for (i in 1:ncol(data.all)) { scatter.smooth(data.all[,i], resid(cox.molec_preop.null, type = "martingale"), main = colnames(data.all)[i]); abline(h = 0) }
# Looks like a threshold effect for size.  Introduce a new
# variable to code this.
data.all$Path.Size.Cent.Thresh = pmin(data.all$Path.Size.Cent, 0)
dd <<- datadist(data.all)
options(datadist = "dd")

# zph analysis indicated violations for Patient.Sex and Molec.S100A4.DCThresh.
# Stratify on these to get around the problem
cox.molec_preop = cph(
	Surv(time, event) ~ 
		strat(Patient.Sex) + 
		History.Diagnosis.AgeAt.Cent + 
		Path.LocationBody + 
		Path.Size.Cent.Thresh + 
		strat(Molec.S100A4.DCThresh) + 
		Molec.S100A2.DCThresh, 
	data = data.all, x = TRUE, y = TRUE, surv = FALSE)
# cox.zph(cox.molec_preop)			# All good now
# par(ask = TRUE)
# for (i in 1:ncol(data.all)) { scatter.smooth(data.all[,i], resid(cox.molec_preop, type = "deviance"), main = colnames(data.all)[i]); abline(h = 0) }

weibull.molec_preop = psm(
	Surv(time, event) ~ 
		Patient.Sex + 
		History.Diagnosis.AgeAt.Cent + 
		Path.LocationBody + 
		Path.Size.Cent + 
		Molec.S100A4.DCThresh + 
		Molec.S100A2.DCThresh, 
	data = data.all, dist = "weibull", x = TRUE, y = TRUE)
# anova(weibull.molec_preop)
# fastbw(weibull.molec_preo)p
# validate(weibull.molec_preop, B = 1000)


# temp.start = 1
# for (temp.i in 1:length(km.molec_preop$strata))
# {
# 	temp.times = km.molec_preop$time[temp.start:(temp.start+km.molec_preop$strata[temp.i])]
# 	temp.surv = km.molec_preop$surv[temp.start:(temp.start+km.molec_preop$strata[temp.i])]
# 	if (temp.i == 1) { 
# 		plot(log(temp.times), log(-log(temp.surv)), type = "s", col = temp.i, xlab = "log(t)", ylab = "log(-log(hat(S)(t)))", ylim = c(-5, 2)) 
# 	} else {
# 		lines(log(temp.times), log(-log(temp.surv)), type = "s", col = temp.i) 
# 	}
# 	temp.start = temp.start + km.molec_preop$strata[temp.i] + 1
# }
# Not great fits to Weibull.  :/


library(RColorBrewer)
pal = brewer.pal(3, "Dark2")


calcSurvivalCurve = function(input, which)
{
	if (which == 1)
	{
		newdata = data.frame(
			Patient.Sex = factor(input$Patient.Sex1, levels = c("M", "F")), 
			History.Diagnosis.AgeAt.Cent = input$History.Diagnosis.AgeAt1 - 68,
			Path.LocationBody = !input$Path.LocationHead1, 
			Path.Size.Cent = input$Path.Size1 - 30, 
			Path.Size.Cent.Thresh = pmin(input$Path.Size1 - 30, 0),
			Molec.S100A2.DCThresh = input$Molec.S100A2.DCThresh1, 
			Molec.S100A4.DCThresh = input$Molec.S100A4.DCThresh1)
		prognostic = input$Prognostic1
	}
	else
	{
		newdata = data.frame(
			Patient.Sex = factor(input$Patient.Sex2, levels = c("M", "F")), 
			History.Diagnosis.AgeAt.Cent = input$History.Diagnosis.AgeAt2 - 68,
			Path.LocationBody = !input$Path.LocationHead2, 
			Path.Size.Cent = input$Path.Size2 - 30, 
			Path.Size.Cent.Thresh = pmin(input$Path.Size2 - 30, 0),
			Molec.S100A2.DCThresh = input$Molec.S100A2.DCThresh2, 
			Molec.S100A4.DCThresh = input$Molec.S100A4.DCThresh2)		
		prognostic = input$Prognostic2
	}

	if (prognostic == "rf")
	{
		temp = predict(rsf.molec_preop, newdata = newdata)
		predobj = list(time = temp$time.interest, surv = temp$survival, lower = rep(NA, length(temp$time.interest)), upper = rep(NA, length(temp$time.interest)))
	}
	else if (prognostic %in% c("cox", "weibull"))
	{
		temp = survest(switch(prognostic, "cox" = cox.molec_preop, "weibull" = weibull.molec_preop), newdata = newdata, conf.int = 0.95)
		predobj = temp
	}
	else if (prognostic == "km")
	{
		temp.stratum = sprintf("Molec.S100A4.DCThresh=%s,Molec.S100A2.DCThresh=%s", newdata$Molec.S100A4.DCThresh, newdata$Molec.S100A2.DCThresh)
		temp.sel = km.molec_preop.strata == temp.stratum
		predobj = list(time = km.molec_preop$time[temp.sel], surv = km.molec_preop$surv[temp.sel], lower = km.molec_preop$lower[temp.sel], upper = km.molec_preop$upper[temp.sel])
	}

	predobj
}


alphaCol = function(col, alpha)
{
	do.call(rgb, c(as.list(col2rgb(col)/255), alpha = alpha))
}


plotSurvCI = function(times, lower, upper, col)
{
	# lines(times, lower, type = "s", lwd = 1, col = col, lty = "dotted")
	# lines(times, upper, type = "s", lwd = 1, col = col, lty = "dotted")
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
