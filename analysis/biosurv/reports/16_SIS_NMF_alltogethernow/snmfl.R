library(nnls)

snmfl.atom = function(A, k, eta, beta, maxiter, converr)
{
	if (eta < 0)	{ eta = max(A) }
	n = nrow(A)
	p = ncol(A)

	W = matrix(runif(n*k)*max(A), nrow = n)
	H = matrix(0, nrow = k, ncol = p)

	P = rbind(W, sqrt(eta) * diag(k))
	Q = rbind(t(H), sqrt(beta) * rep(1, k))

	lasterr = Inf
	conv = FALSE
	for (i in 1:maxiter)
	{
		thisXhat = W %*% H
		thisdiff = thisXhat - A
		thiserr = sum(thisdiff[,]^2)

		if (lasterr - thiserr < converr)
		{
			conv = TRUE
			break
		}
		lasterr = thiserr

		P[1:nrow(W),] = W
		for (j in 1:p)
		{
			H[,j] = coefficients(nnls(P, c(A[,j], rep(0, k))))
		}

		Q[1:ncol(H),] = t(H)
		for (j in 1:n)
		{
			W[j,] = t(coefficients(nnls(Q, c(A[j,], 0))))
		}
	}

	list(W = W, H = H, i = i, converged = conv)
}


snmfl = function(A, ks, nrun = 1, eta = -1, beta = 0.01, maxiter = 1e3, converr = 1e-3, cores = 1)
{
	jobs = expand.grid(k = ks, run = 1:nrun)

	if (cores == 1)
	{
		results = lapply(1:nrow(jobs), function(job_i) { cat(job_i, "/", nrow(jobs), "\n"); snmfl.atom(A, jobs$k[job_i], eta, beta, maxiter, converr) } )
	}
	else
	{
		require(parallel)
		results = mclapply(1:nrow(jobs), function(job_i) { cat(job_i, "/", nrow(jobs), "\n"); snmfl.atom(A, jobs$k[job_i], eta, beta, maxiter, converr) }, mc.cores = cores)
	}

	results_for_each_k = tapply(1:nrow(jobs), jobs$k, function(job_is) {
		this_k_results = results[job_is]
		norms = sapply(this_k_results, function(r) sum((A - (r$W %*% r$H))^2))
		best_norm = which.min(norms)
		list(best_fit = this_k_results[[best_norm]], norms = norms)
	})
	results_for_each_k = results_for_each_k[order(as.numeric(names(results_for_each_k)))]

	results_for_each_k
}


# basis == W
# coef == H
