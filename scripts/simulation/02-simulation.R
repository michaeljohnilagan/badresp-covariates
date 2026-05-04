# import packages
library('R6')
library('MASS')
library('detranli')

# import my stuff
ao = new.env()
with(ao, {
	source('./utils/aomodels.R', local=TRUE)
	source('./utils/classification.R', local=TRUE)
})

# initialize objects
settings = list()
samplers = list()

# settings: L1P1
settings$l1p1 = with(new.env(), {
	feat_funs = c('mahal', 'ptcossim') # nonresponsivity indices
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
	list(cbt=cbt)
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
	var_noise = 1 # noise variance in each feature
	numcols_total = 3 # how many covariates, good + bad
	list(var_noise=var_noise, numcols_total=numcols_total)
})

# samplers: covariates
Sys.time(); samplers$covariates = with(new.env(), {
	sampler = function(labels, numcols_useful) {
		stopifnot(numcols_useful<=settings$covariates$numcols_total)
		numcols_useless = settings$covariates$numcols_total-
		numcols_useful
		loadings = c(rep(1, times=numcols_useful), 
		rep(0, times=numcols_useless))
		systematic = outer(labels, loadings)
		noise_num_cells = length(labels)*
		settings$covariates$numcols_total
		noise = rnorm(noise_num_cells, 0, 
		sqrt(settings$covariates$var_noise))
		systematic+noise
	}
	# output
	list(sampler=sampler)
})

# function: run simulation replicate, sampler
run_repl_sampler = function(n, contam, xnumuf, zwhich) {
	# true class label
	y = ifelse((1:n)<=(n*contam), 1, 0)
	# generate likert
	included = settings$likert$inclusion[[zwhich]]
	z = t(sapply(1:n, function(i) {
		if(y[i]==1) {
			samplers$likert[['cnr']]()
		} else if(y[i]==0) {
			samplers$likert[['noncnr']]()
		} else {
			stop('invalid class label')
		}
	}))[, included]
	# generate covariates
	x = samplers$covariates$sampler(labels=y, numcols_useful=xnumuf)
	# put together
	return(list(y=y, z=z, x=x))
}

# function: run simulation replicate
run_repl = function(n, contam, xnumuf, zwhich, verbose=FALSE) {
	# generate data
	dat = run_repl_sampler(n=n, contam=contam, xnumuf=xnumuf, 
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
	pred_class_lab_ao0 = round(mod0$calc_postr()[, '1'])
	# predictions: AO1 classifier
	init = c(mod0$par_get()$steepness, qlogis(mod0$par_get()$prevalence), 
	rep(0, times=num_features-1))
	mod1 = ao$AO1Model$new(tables=settings$ao$cbt)
	trymod1 = try(mod1$fit(success_counts=succ, features=x, 
	init=init), silent=TRUE) # canned ML
	# default to AO0 if AO1 fails
	optim_exception = class(trymod1)=='try-error'
	if(optim_exception) {
		mod1 = mod0$clone()
	}
	pred_class_lab_ao1 = round(mod1$calc_postr()[, '1'])
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
if(TRUE) {
	set.seed(91)
	foo = run_repl(n=200, contam=0.75, xnumuf=1, zwhich='even', 
	verbose=TRUE)
	print(with(foo, {
		do.call(rbind, list(sc=met_sc, ao0=met_ao0, ao1=met_ao1))
	}), digits=3)
	rm(foo)
}

# function: run simulation cell
run_cell = function(numrepl, n, contam, xnumuf, zwhich) {
	# scenario information
	scenario = data.frame(n=n, contam=contam, xnumuf=xnumuf, 
	zwhich=zwhich)
	# do many replicates
	start_seed = 1000
	result_cell = lapply(start_seed+1:numrepl, function(seed) {
		set.seed(seed)
		result_repl = do.call(run_repl, scenario)
		unlist(result_repl[c('met_sc', 'met_ao0', 'met_ao1', 
		'll_improved', 'optim_exception')])
	})
	# put together
	simplified = as.data.frame(do.call(rbind, result_cell))
	together = cbind(scenario, simplified)
	return(together)
}

# demonstrate one cell
if(FALSE) {
	set.seed(99)
	Sys.time(); foo = run_cell(numrepl=3, n=200, contam=0.75, 
	xnumuf=1, zwhich='even'); Sys.time()
	print(foo, digits=3)
	rm(foo)
}

# simulation study factors
sim_factors = list(n=c(100, 300, 900), contam=c(5, 25, 50, 75, 95)/100,
xnumuf=0:3, zwhich=c('even', 'all'))
print(sim_factors)

# simulation number of replicates
numrepl = 3
message('this run has ', numrepl, ' replicates per cell')

# simulation results
sim_results = array(list(), dim=sapply(sim_factors, length))
for(i1 in 1:length(sim_factors$n))
for(i2 in 1:length(sim_factors$contam))
for(i3 in 1:length(sim_factors$xnumuf))
for(i4 in 1:length(sim_factors$zwhich)) {
	with(sim_factors, {
		scenario_id = mapply('[[', sim_factors, c(i1, i2, i3, i4))
		message(format(Sys.time(), "%Y-%m-%d %H:%M"), ' | ', 
		paste(scenario_id, collapse=' '))
	}) # report cell
	sim_results[[i1, i2, i3, i4]] = with(sim_factors, {
		run_cell(numrepl=numrepl, n=n[i1], contam=contam[i2], 
		xnumuf=xnumuf[i3], zwhich=zwhich[i4])
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

