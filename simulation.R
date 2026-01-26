# import packages
library('R6')
library('MASS')
library('detranli')

# import my stuff
Sys.time(); ao = new.env()
with(ao, {
	source('./utils/aomodels.R', local=TRUE)
	source('./utils/classification.R', local=TRUE)
}); Sys.time()

# initialize objects
settings = list()
samplers = list()

# settings: L1P1
settings$l1p1 = with(new.env(), {
	feat_funs = c('mahal', 'ptcor') # nonresponsivity indices
	feat_idvals = c(0, +1) # ideal point
	numperms = 200 # number of synthetic rows per observed
	nomsens = 0.95 # nominal sensitivity
	list(feat_funs=feat_funs, feat_idvals=feat_idvals, numperms=numperms, 
	nomsens=nomsens)
})

# settings: AO models
Sys.time(); settings$ao = with(new.env(), {
	numperms = settings$l1p1$numperms
	shift_limit = -10 # most negative shift possible to fit
	num_gridpoints = 300 # number of grid points to compute
	cbt = ao$CarpBinTable$new(size=numperms, 
	shift_limit=shift_limit, num_gridpoints=num_gridpoints) # table set
	em_maxit = 50 # maximum iterations for EM algorithm
	list(cbt=cbt, em_maxit=em_maxit)
}); Sys.time()

# settings: likert
settings$likert = with(new.env(), {
	# known facts about the dataset
	filename = './data/openpsychometrics-hsq.csv'
	column_idx_likert = 1:32
	code_missing = -1
	pointscales = rep(5, times=length(column_idx_likert))
	# load dataset and process
	dat = read.csv(filename, header=TRUE)[, column_idx_likert]
	dat = replace(dat, dat==code_missing, NA) # code missing
	rows_invalid_bool = apply(dat, 1, function(r) {
		all(is.na(r))
	}) # rows to be removed
	dat = as.matrix(subset(dat, !rows_invalid_bool)) # remove empty rows
	# inclusion definitions for zwhich
	numcols = ncol(dat)
	inclusion = list(all=rep(TRUE, times=numcols), 
	even=((1:numcols)%%2)==0)
	# output
	list(dat=dat, pointscales=pointscales, inclusion=inclusion)
})

# samplers: likert
samplers$likert = with(new.env(), {
	sampler_noncnr = function() {
		detranli::samplerows(1, data=settings$likert$dat)
	}
	sampler_cnr = function() {
		detranli::rcnrbinom(1, pointscales=settings$likert$pointscales)
	}
	list(noncnr=sampler_noncnr, cnr=sampler_cnr)
})

# settings: covariates
settings$covariates = with(new.env(), {
	usefulcorr = 0 # correlation between useful covariates
	numcols_useful = 2 # how many useful covariates
	numcols_useless = 8 # how many useless covariates
	xsep_auc = setNames(c(0.5, 0.6, 0.7, 0.8), 
	c('none', 'small', 'medium', 'large')) # target AUCs
	list(usefulcorr=usefulcorr, numcols_useful=numcols_useful, 
	numcols_useless=numcols_useless, xsep_auc=xsep_auc)
})

# samplers: covariates
Sys.time(); samplers$covariates = with(new.env(), {
	# use the settings
	usefulcorr = settings$covariates$usefulcorr
	numcols_useful = settings$covariates$numcols_useful
	numcols_useless = settings$covariates$numcols_useless
	numcols = numcols_useful+numcols_useless
	targets = settings$covariates$xsep_auc
	# class 0 mean and common variance
	mu0 = rep(0, times=numcols)
	usefulness = ifelse((1:numcols)<=numcols_useful, TRUE, FALSE)
	sig = matrix(NA, numcols, numcols)
	for(i in 1:numcols) for(j in 1:numcols) {
		if(i==j) {
			sig[i, j] = 1
		} else if(usefulness[i]==usefulness[j]) {
			sig[i, j] = usefulcorr
		} else {
			sig[i, j] = 0
		}
	}
	# AUC as a function of the shift
	shift2auc = function(shift) {
		# implied class 1 mean
		mu1 = c(rep(shift, times=numcols_useful), rep(shift, 
		times=numcols_useless))
		# implied noncentrality parameter
		ncp = as.vector(t(mu1)%*%solve(sig)%*%mu1)
		# integration for AUC
		spec2sens = function(spec) {
			qtile = qchisq(spec, df=numcols, ncp=0)
			1-pchisq(qtile, df=numcols, ncp=ncp)
		}
		integrate(spec2sens, lower=0, upper=0.999)$value
	}
	# shift as a function of the AUC
	auc2shift = function(auc) {
		optim(0.5, function(s) {
			(shift2auc(s)-auc)^2
		}, method='L-BFGS-B', lower=0, upper=2)$par
	}
	# for each target AUC, get the shift
	shifts = setNames(c(0, sapply(targets[-1], auc2shift)), 
	names(targets))
	# class 1 mean by target
	mu1 = lapply(shifts, function(s) {
		c(rep(s, times=numcols_useful), rep(0, 
		times=numcols_useless))
	})
	# print the implied parameters
	print(list(mu0=mu0, mu1=mu1, sig=sig))
	# generators
	sampler_noncnr = function() {
		MASS::mvrnorm(1, mu=mu0, Sig=sig)
	}
	sampler_cnr = function(xsep) {
		MASS::mvrnorm(1, mu=mu1[[xsep]], Sig=sig)
	}
	# output
	list(noncnr=sampler_noncnr, cnr=sampler_cnr)
}); Sys.time()

# demonstrate xsep in 2D
set.seed(840)
pdf('xsep.pdf')
invisible(sapply(names(settings$covariates$xsep_auc), function(u) {
	# true class label
	y = rep(0:1, times=1e3)
	# coordinates
	x = t(sapply(y, function(yy) {
		if(yy==1) {
			samplers$covariates[['cnr']](u)
		} else if(yy==0) {
			samplers$covariates[['noncnr']]()
		}
	}))
	# colors
	opacity = 0.2
	col0 = rgb(0, 0, 1, alpha=opacity)
	col1 = rgb(1, 0, 0, alpha=opacity)
	col = ifelse(y==1, col1, col0)
	# plot
	plot(x[, 1], x[, 2], col=col, pch=19, cex=1.5,
	xlab='covariate 1', ylab='covariate 2', main=u)
}))
dev.off()

# function: run simulation replicate, sampler
run_repl_sampler = function(n, contam, xsep, zwhich) {
	# true class label
	y = ifelse((1:n)<(n*contam), 1, 0)
	# generate likert
	included = settings$likert$inclusion[[zwhich]]
	z = t(sapply(y, function(yy) {
		if(yy==1) {
			samplers$likert[['cnr']]()
		} else if(yy==0) {
			samplers$likert[['noncnr']]()
		}
	}))[, included]
	# generate covariates
	x = t(sapply(y, function(yy) {
		if(yy==1) {
			samplers$covariates[['cnr']](xsep)
		} else if(yy==0) {
			samplers$covariates[['noncnr']]()
		}
	}))
	# put together
	return(list(y=y, z=z, x=x))
}

# function: run simulation replicate
run_repl = function(n, contam, xsep, zwhich, verbose=FALSE) {
	# generate data
	dat = run_repl_sampler(n=n, contam=contam, xsep=xsep, 
	zwhich=zwhich)
	included = settings$likert$inclusion[[zwhich]]
	# use settings for L1P1
	pointscales = settings$likert$pointscales[included]
	numperms = settings$l1p1$numperms
	feat_funs = settings$l1p1$feat_funs
	feat_idvals = settings$l1p1$feat_idvals
	nomsens = settings$l1p1$nomsens
	# calculate p values
	p = detranli::cnrdetect(dat$z, pointscales=pointscales, 
	numperms=numperms, feat_funs=feat_funs, feat_idvals=feat_idvals)
	p = ifelse(is.na(p), 1, p)
	succ = ao$pval2count(p, size=numperms)
	# add intercept to features
	x = cbind(1, dat$x)
	num_features = ncol(x)
	# predictions: sensitivity calibrated (SC) classifier
	pred_class_lab_sc = ifelse(p>=(1-nomsens), 1, 0)
	# predictions: AO0 clasiffier
	mod0 = ao$AO0Model$new(tables=settings$ao$cbt)
	mod0$fit_mm(success_counts=succ)
	pred_class_lab_ao0 = round(mod0$calc_postr_cnr())
	# predictions: AO1 classifier
	init = c(mod0$par_get()$steepness, qlogis(mod0$par_get()$prevalence), 
	rep(0, times=num_features-1))
	mod1 = ao$AO1Model$new(tables=settings$ao$cbt)
	if(FALSE) {
		trymod1 = try(mod1$fit_em(success_counts=succ, features=x, 
		init=init,maxit=settings$ao$em_maxit, 
		verbose=verbose), silent=TRUE) # EM algorithm
	} else {
		trymod1 = try(mod1$fit(success_counts=succ, features=x, 
		init=init), silent=TRUE) # canned ML
	}
	# default to AO0 if AO1 fails
	optim_exception = class(trymod1)=='try-error'
	if(optim_exception) {
		mod1 = mod0$clone()
	}
	pred_class_lab_ao1 = round(mod1$calc_postr_cnr())
	# get metrics
	metrics = function(predicted_class_label) {
		ao$metrics(true_class_label=dat$y, 
		predicted_class_label=predicted_class_label)
	}
	met = lapply(list(met_sc=pred_class_lab_sc, 
	met_ao0=pred_class_lab_ao0, met_ao1=pred_class_lab_ao1), metrics)
	# append log likelihood to metrics
	ll_sc = setNames(NA, 'll')
	ll_ao0 = setNames(mod0$loglikelihood(), 'll')
	ll_ao1 = setNames(mod1$loglikelihood(), 'll')
	ll_improved = unname(sign(ll_ao1-ll_ao0))
	met = mapply(c, met, list(ll_sc, ll_ao0, ll_ao1), SIMPLIFY=FALSE)
	# put together
	together = c(list(y=dat$y, z=dat$z, x=x, p=p, succ=succ, 
	ll_improved=ll_improved, optim_exception=optim_exception), met)
	return(together)
}

# demonstrate single replicate
set.seed(91)
foo = run_repl(n=200, contam=0.75, xsep='small', zwhich='even', 
verbose=TRUE)
print(with(foo, {
	do.call(rbind, list(sc=met_sc, ao0=met_ao0, ao1=met_ao1))
}), digits=3)
rm(foo)

# function: run simulation cell
run_cell = function(numrepl, n, contam, xsep, zwhich) {
	# scenario information
	scenario = data.frame(n=n, contam=contam, xsep=xsep, 
	zwhich=zwhich)
	# do many replicates
	start_seed = 1000
	result_cell = lapply(start_seed+1:numrepl, function(seed) {
		set.seed(seed)
		result_repl = run_repl(n=n, contam=contam, 
		xsep=xsep, zwhich=zwhich)
		c(unlist(result_repl[c('met_sc', 'met_ao0', 'met_ao1', 
		'll_improved', 'optim_exception')]))
	})
	# put together
	simplified = as.data.frame(do.call(rbind, result_cell))
	together = cbind(scenario, simplified)
	return(together)
}

# demonstrate one cell
set.seed(99)
Sys.time(); foo = run_cell(numrepl=3, n=200, contam=0.75, xsep='small', 
zwhich='even'); Sys.time()
print(foo, digits=3)
rm(foo)

# simulation study factors
sim_factors = list(n=c(100, 300, 900), contam=c(5, 25, 50, 75, 95)/100,
xsep=c('none', 'small', 'medium', 'large'), 
zwhich=c('even', 'all'))
print(sim_factors)

# simulation number of replicates
numrepl = 3
print(numrepl)

# simulation results
sim_results = array(list(), dim=sapply(sim_factors, length))
for(i1 in 1:length(sim_factors$n))
for(i2 in 1:length(sim_factors$contam))
for(i3 in 1:length(sim_factors$xsep))
for(i4 in 1:length(sim_factors$zwhich)) {
	with(sim_factors, {
		scenario_id = mapply('[[', sim_factors, c(i1, i2, i3, i4))
		message(format(Sys.time(), "%Y-%m-%d %H:%M"), ' | ', 
		paste(scenario_id, collapse=' '))
	}) # report cell
	sim_results[[i1, i2, i3, i4]] = with(sim_factors, {
		run_cell(numrepl=numrepl, n=n[i1], contam=contam[i2], 
		xsep=xsep[i3], zwhich=zwhich[i4])
	}) # save cell result
}; Sys.time()
rm(i1, i2, i3, i4)
warns = warnings() # save warnings
print(warns)

# function: summarize cell
summarize_cell = function(tab) {
	# get scenario
	param_names = names(sim_factors)
	scenario = tab[1, param_names]
	# get outcome measures summary
	metric_indices = setdiff(colnames(tab), param_names)
	outcome_measures = colMeans(tab[, metric_indices])
	return(as.data.frame(c(scenario, outcome_measures)))
}

# summarize cells
sim_tab = do.call(rbind, lapply(sim_results, summarize_cell))

# end session
Sys.time()
save.image("./simulation.RData")
devtools::session_info()
