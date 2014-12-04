breakTies = function(s)
{
	if (!any(duplicated(s[,1])))
	{
		return(s)
	}

	deltas = diff(sort(s[,1]))
	mindelta = min(deltas[zapsmall(deltas) > 0])

	newtimes = s[,1]
	while (any(duplicated(newtimes)))
	{
		newtimes = s[,1] + (runif(length(newtimes))-0.5)*mindelta*duplicated(s[,1])
	}

	Surv(newtimes, s[,2])
}


SIS.FAST = function(x, y, gamma, scale = FALSE, ...)
{
	require(ahaz)
	if (scale == TRUE)
	{
		xstd = (x - rowMeans(x)) / (apply(x, 1, sd))
	}
	else
	{
		xstd = x
	}
	fast_scores = ahaz(breakTies(y), t(xstd), univariate = TRUE)$d
	abs(fast_scores) > gamma
}


# autoThresh = function(x, y, selfunc, Js, B, ...)
# {
# 	message("Searching for best selection threshold...")
# 	p = nrow(x)
# 	scores = sapply(Js, function(j) {
# 		message(sprintf("  Testing %d", j))
# 		obs = selfunc(x, y, n = j, ...)
# 		boots = sapply(1:B, function(bi) {
# 			samp = sample.int(ncol(x), replace = TRUE)
# 			selfunc(x[,samp], y[samp,], n = j, ...)
# 		})
# 		overlap = obs & boots
# 		noverlap = colSums(overlap)
# 		EH0 = (j^2)/p
# 		varH0 = (j^2*(p-j)^2)/(p^2*(p-1))

# 		score = (mean(noverlap) - EH0) / sqrt(varH0/B)
# 		score
# 	})
# 	best_score = which.max(scores)
# 	plot(scores ~ Js, xlab = "Number of variables", ylab = "Selection score", type = "o")
# 	abline(v = Js[best_score])
# 	message(sprintf("Best score found for j = %d", Js[best_score]))
# 	selfunc(x, y, n = Js[best_score])
# }


CPSS = function(x, y, selfunc, tau, B = 50, ...)
{
	n = ncol(x)

	samps = lapply(1:B, function(i) sort(sample.int(n, floor(n/2))))

	Shats = sapply(1:(2*B), function(i) {
		A = samps[[ceiling(i/2)]]
		if (i %% 2 == 0)	{ A = setdiff(1:n, A) }

		selfunc(x[,A], y[A,], ...)
	})

	Pihat = rowSums(Shats) / (2*B)
	qs = colSums(Shats)

	result = list(sel = Pihat > tau, Pihat = Pihat, samples = samps, qs = qs, qhat = mean(qs))
	result
}


corPlot = function(cors, zlim = c(-1, 1), pal = "RdYlGn", scores = NULL, ...)
{
	clust = hclust(dist(cors))
	thepal = brewer.pal(brewer.pal.info[pal,]$maxcolors, pal)

	z = cors[rev(clust$order), clust$order]
	if (!is.null(scores)) {
		scores = t(scores)
		scores = (scores - apply(scores, 1, min)) / as.vector(diff(apply(scores, 1, range)))
		scores = t(apply(scores, 1, function(x) { if(mean(x) < 0.5) x else 1-x } ))
		scores = scores * (zlim[2] - zlim[1]) + zlim[1]
		scores = t(scores)
		scores = scores[,ncol(scores):1]
		z = cbind(z, scores[rev(clust$order),])
	}

	pars = par(no.readonly = TRUE)
	par(mar = c(6, 3, 5, 3)/1.5)
	layout(matrix(c(1, 2), nrow = 1), widths = c(8, 1))
	image(z = z, zlim = zlim, col = thepal, xaxt = "n", yaxt = "n", ...)
	par(mar = c(6, 2, 5, 1)/1.5)
	image(x = c(0, 1), y = seq(zlim[1], zlim[2], length.out = 100), z = matrix(seq(zlim[1], zlim[2], length.out = 99), nrow = 1), col = thepal, xaxt = "n", xlab = "", ylab = "", useRaster = TRUE)
	par(pars)
}

