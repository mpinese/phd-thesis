load("01_NSWPCN.rda")

#library(gplots)
#heatmap.2(1 - 1*is.na(data), trace = "none", useRaster = TRUE)

library(survival)

fit = survfit(Surv(History.Death.EventTimeDays, History.DSDeath.Event) ~ Treat.Resected, data)
plot(fit, xlim = c(0, 3000), col = 2:1)
data$Molec.S100A4.DCThresh = (data$Molec.S100A4.NucInt >= "1") | (data$Molec.S100A4.CytoInt >= "1")
data$Molec.S100A2.DCThresh = (data$Molec.S100A2.Int >= "3") & (data$Molec.S100A2.Percent > 30)
data$group = ifelse(data$Treat.Resected, data$Molec.S100A4.DCThresh*2 + data$Molec.S100A2.DCThresh + 1, 0)
fit = survfit(Surv(History.Death.EventTimeDays, History.DSDeath.Event) ~ group, data)
plot(fit, xlim = c(0, 3000), col = 2:6)
legend("topright", legend = c("Biopsy", "A2- A4-", "A2+ A4-", "A2- A4+", "A2+ A4+"), col = 2:6, lty = "solid")


approxHaz = function(fit, df)
{
	nstrata = length(fit$strata)
	strata = rep(names(fit$strata), fit$strata)
	strata.times = tapply(1:length(strata), strata, function(is) fit$time[is])
	strata.surv = tapply(1:length(strata), strata, function(is) fit$surv[is])
	strata.cumhaz = lapply(strata.surv, function(x) -log(x))

	strata.cumhaz.splines = mapply(function(times, cumhaz) {
		sel = times != 0 & cumhaz != 0
		times = times[sel]
		cumhaz = cumhaz[sel]
		smooth.spline(log(times), log(cumhaz), df = df)
	}, strata.times, strata.cumhaz, SIMPLIFY = FALSE)

	delta = 0.1
	grid = seq(0, max(fit$time), delta)
	strata.haz_estimate = lapply(strata.cumhaz.splines, function(spl) {
		cumhaz = exp(predict(spl, x = log(grid))$y)
		haz = diff(c(0, cumhaz)) / delta
		haz = pmax(0, haz)
		cumhaz_corrected = cumsum(haz * delta)
		#surv_corrected = exp(-cumhaz_corrected)
		haz_corrected = diff(c(0, cumhaz_corrected)) / delta
		haz_corrected
	})

	plot(0 ~ 0, type = "n", xlim = c(1, max(fit$time)), ylim = c(0, 1))
	for (i in 1:nstrata)
	{
		lines(strata.surv[[i]] ~ strata.times[[i]], type = "s", col = i, lwd = 2)
		temp.cumhaz_corrected = cumsum(strata.haz_estimate[[i]] * delta)
		temp.surv_corrected = exp(-temp.cumhaz_corrected)
		lines(grid, temp.surv_corrected, lwd = 1, col = i, lty = "dashed")
	}

	plot(0 ~ 0, type = "n", xlim = c(1, 2000), ylim = c(0, max(sapply(strata.haz_estimate, max))))
	for (i in 1:nstrata)
	{
		lines(grid, strata.haz_estimate[[i]], lwd = 2, col = i, lty = "solid")
	}
	legend("topright", legend = paste(names(fit$strata), fit$strata, sep = " "), col = 1:nstrata, lty = "solid", lwd = 2)

	strata.haz_estimate
}


temp.y = Surv(data$History.Death.EventTimeDays, data$History.DSDeath.Event)
temp.sel = !is.na(temp.y[,1]) & !is.na(temp.y[,2]) & data$Treat.Resected == TRUE

temp.y = temp.y[temp.sel,]
temp.x = data[temp.sel, !(colnames(data) %in% c("History.Death.EventTimeDays", "History.DSDeath.Event", "Patient.ID"))]
temp.x$Molec.S100A4.CytoH = temp.x$Molec.S100A4.CytoPercent * as.numeric(as.character(temp.x$Molec.S100A4.CytoInt))
temp.x$Molec.S100A4.NucH = temp.x$Molec.S100A4.NucPercent * as.numeric(as.character(temp.x$Molec.S100A4.NucInt))
temp.x$Molec.Grb7.H = temp.x$Molec.Grb7.Percent * as.numeric(as.character(temp.x$Molec.Grb7.Int))
temp.x$Molec.HENT1.H = temp.x$Molec.HENT1.Percent * as.numeric(as.character(temp.x$Molec.HENT1.Int))
temp.x$Molec.HOXB2.H = temp.x$Molec.HOXB2.Percent * as.numeric(as.character(temp.x$Molec.HOXB2.Int))
temp.x$Molec.S100A2.H = temp.x$Molec.S100A2.Percent * as.numeric(as.character(temp.x$Molec.S100A2.Int))

temp.x$Molec.S100A4.DCThresh = (temp.x$Molec.S100A4.NucInt >= "1") | (temp.x$Molec.S100A4.CytoInt >= "1")
temp.x$Molec.S100A2.DCThresh = (temp.x$Molec.S100A2.Int >= "3") & (temp.x$Molec.S100A2.Percent > 30)

surv.pvals = sapply(temp.x, function(x1) { cat("."); res = try(summary(coxph(temp.y ~ x1))$logtest[["pvalue"]]); if (class(res) == "try-error") { return(NA) } else { return(res) } } )
sort(surv.pvals)

# Dates of death, diagnosis, and recurrence, are highly significant.  Possibilities:
#   * Change in survival over time due to changes in management
#     -- this could be explained by values of other fields (eg chemo changes)
#     -- or could be other unmeasured variables (eg surgical management)
#   * Change in cohort composition over time, possibly due to different recruitment
#     priorities / centres.
#   * A bit of both.

# S100A2, A4 etc are significant, but only when Path.Subtype contains all values.
# Restricting to NotSpecified removes these predictors.  Could it be that the 
# molecular measures are in fact surrogate indicators of subtype, which is what
# is actually determining outcome?  This will be hard to check with these data,
# as the majority of samples have a subtype of NotSpecified -- presumably many
# of these are PDAC, but I can't be sure of this.

subtype.pvals = sapply(temp.x[,colnames(temp.x) != "Path.Subtype"], function(x1) {
	if ("factor" %in% class(x1) || "logical" %in% class(x1))
	{
		if (length(unique(x1)) <= 1)	{ return(NA) }
		print(table(x1))
		pval = try(fisher.test(temp.x$Path.Subtype, x1)$p.value)
		if (class(pval) == "try-error")
		{
			pval = chisq.test(temp.x$Path.Subtype, x1)$p.value
		}
		return(pval)
	}
	else if ("numeric" %in% class(x1))
	{
		return(anova(lm(x1 ~ temp.x$Path.Subtype))["Pr(>F)"][1,])
	}
	else
	{
		return(NA)
	} })
sort(subtype.pvals)

boxplot(Molec.S100A2.H ~ Path.Subtype, temp.x)
boxplot(Molec.S100A2.StromaScore ~ Path.Subtype, temp.x)
plot(table(temp.x$Molec.S100A2.DCThresh, temp.x$Path.Subtype))
plot(table(temp.x$Molec.S100A4.DCThresh, temp.x$Path.Subtype))
plot(table(temp.x$Treat.Chemo.Palliative, temp.x$Path.Subtype))

# Yerrrp.  I think the markers are definitely reflecting subtype.
# Chemo makes sense, seeing as it'd be a postop decision based on path.

# So there's a bit of a problem.  The results per-se are fine, because
# this nomogram is intended to cover all-comers pre-surgery, before the 
# particular subtype would ordinarily be known.  The idea is to take
# a biopsy of the Ca prior to resection, assess the markers, and then 
# determine likely post-operative survival based on that.  But if most
# of the information is just encoded in subtype anyway... why bother?
# Originally the biopsies were to be FNAs, where insufficient undamaged 
# tissue would be present to make a call on subtype.  But technically
# it's much simpler to stain and score cores, and practically I think this
# is where the technique is heading.  In that case, could subtype be
# effectively inferred from core anyway?

# Read DC's paper to see if he covers this.

# DC's paper hinges on EUS-FNA.  So subtype is out.  But do mention
# the link to DC; he may find it interesting.

# DC's paper looks at A2 and A4 only.  Specifically, he used the following
# derived variables:
#   A2: binary.  Cytoplasmic staining of 3+ in > 30% of cells.
#   A4: binary.  Cytoplasmic OR nuclear staining of any intensity in > 1% of cells.

# There may be some leeway to optimize these, especially considering
# the large number of validation cohorts now accessible.

sort(sapply(temp.x, function(x1) sum(!is.na(x1))), decreasing = TRUE)
sort(sapply(temp.x, function(x1) sum(!is.na(x1))), decreasing = TRUE) / nrow(temp.x)

# Subset to see just how many we've got with the required data.
data$Molec.S100A4.DCThresh = (data$Molec.S100A4.NucInt >= "1") | (data$Molec.S100A4.CytoInt >= "1")
data$Molec.S100A2.DCThresh = (data$Molec.S100A2.Int >= "3") & (data$Molec.S100A2.Percent > 30)
temp.sel = 
	!is.na(data$History.Death.EventTimeDays) & 
	!is.na(data$History.DSDeath.Event) & 
	data$Treat.Resected == TRUE & 
	!is.na(data$Molec.S100A2.DCThresh) &
	!is.na(data$Molec.S100A4.DCThresh)
sum(temp.sel)
mean(temp.sel)

# 314 patients (~30%).  Doesn't add up with DC's paper numbers, but oh well.

