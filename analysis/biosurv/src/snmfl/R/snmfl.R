snmfl = function(A, ks, nrun = 1, eta = -1, alpha = 0.01, max_iter = 1e3, conv_delta = 1e-3, conv_interval = 1, cores = 1)
{
	jobs = expand.grid(k = ks, run = 1:nrun)

	if (cores == 1)
	{
		results = lapply(1:nrow(jobs), function(job_i) { cat(job_i, "/", nrow(jobs), "\n"); snmfl_atomic(A, jobs$k[job_i], eta, alpha, max_iter, conv_delta, conv_interval) } )
	}
	else
	{
		require(parallel)
		results = mclapply(1:nrow(jobs), function(job_i) { cat(job_i, "/", nrow(jobs), "\n"); snmfl_atomic(A, jobs$k[job_i], eta, alpha, max_iter, conv_delta, conv_interval) }, mc.cores = cores)
	}

	results_for_each_k = tapply(1:nrow(jobs), jobs$k, function(job_is) {
		this_k_results = results[job_is]
		norms = sapply(this_k_results, function(r) r$norm)
		convflags = sapply(this_k_results, function(r) r$converged)
		best_norm = which.min(norms)
		list(best_fit = this_k_results[[best_norm]], norms = norms, convergence_flags = convflags)
	})
	results_for_each_k = results_for_each_k[order(as.numeric(names(results_for_each_k)))]

	results_for_each_k
}
