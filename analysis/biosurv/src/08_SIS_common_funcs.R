breakTies = function(s)
{
	while (any(duplicated(s[,1])))
	{
		ts = sort(s[,1])
		deltas = diff(ts)
		deltas = deltas[deltas > diff(range(deltas))/1e4]
		mindelta = min(deltas)

		ties = duplicated(s[,1])
		newtimes = s[,1]
		newtimes[ties] = newtimes[ties] + (runif(sum(ties))-0.5)*mindelta
		s = Surv(newtimes, s[,2])
	}

	s
}


ISIS.FAST = function(x, y, ...)
{
	require(ahaz)
	sel = rep(FALSE, nrow(x))
	sel[ahazisis(breakTies(y), t(x), ...)$ISISind] = TRUE
	sel
}


SIS.FAST = function(x, y, ...)
{
	require(ahaz)
	sel = rep(FALSE, nrow(x))
	sel[ahazisis(breakTies(y), t(x), ..., do.isis = FALSE)$SISind] = TRUE
	sel
}


CPSS = function(x, y, selfunc, tau, B = 50, ...)
{
	n = ncol(x)
	p = nrow(x)

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
	par(mar = c(6, 3, 5, 3))
	layout(matrix(c(1, 2), nrow = 1), widths = c(8, 1))
	image(z = z, zlim = zlim, col = thepal, xaxt = "n", yaxt = "n", useRaster = TRUE, ...)
	par(mar = c(6, 2, 5, 1))
	image(x = c(0, 1), y = seq(zlim[1], zlim[2], length.out = 100), z = matrix(seq(zlim[1], zlim[2], length.out = 99), nrow = 1), col = thepal, xaxt = "n", xlab = "", ylab = "", useRaster = TRUE)
	par(pars)
}
