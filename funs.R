# import my stuff
Sys.time()
source('./utils/aomodels.R')
source('./utils/classification.R')
Sys.time()

# test: rcarpal vs dcarpal
set.seed(457)
with(new.env(), {
	# parameters
	sampsize = 4e3
	shift = -1
	# generate and plot
	x = rcarpal(sampsize, shift)
	hist(x, freq=FALSE)
	curve(dcarpal(x, shift=shift), from=0, to=1, lwd=2, add=TRUE)
})

# test: `rcarpbin` vs `dcarpbin`
set.seed(521)
with(new.env(), {
	# parameters
	sampsize = 5e3
	shift = -0.9
	size = 15
	# work generate
	x = rcarpbin(sampsize, size=size, shift=shift)
	masspoints = 0:size
	# calculate PMFs
	pmf_empirical = setNames(sapply(masspoints, function(v) {
		mean(x==v)
	}), masspoints)
	pmf_theoretical = setNames(dcarpbin(masspoints, size=size, shift=shift), 
	masspoints)
	# compare empirical vs theoretical
	plot(masspoints, pmf_empirical, type='h', lwd=2)
	points(masspoints, pmf_theoretical, pch=1, lwd=2)
	round(pmf_empirical-pmf_theoretical, 3)
})

# generate from AO0 hierarchy
rao0 = function(n, size, steepness, prevalence) {
	generated_class0 = rcarpbin(n, size=size, shift=steepness)
	generated_class1 = rcarpbin(n, size=size, shift=0)
	class_label = sample(0:1, size=n, prob=c(1-prevalence, prevalence), 
	replace=TRUE)
	success_count = ifelse(class_label==1, generated_class1, generated_class0)
	return(success_count)
}

# AO0 hierarchy PMF
dao0 = function(x, size, steepness, prevalence) {
	pmf_class0 = dcarpbin(x, size=size, shift=steepness)
	pmf_class1 = dcarpbin(x, size=size, shift=0)
	pmf_mix = prevalence*pmf_class1+(1-prevalence)*pmf_class0
	return(pmf_mix)
}

# test: `dao0` vs `dcarpbin`
with(new.env(), {
	# parameters
	shift = -1
	size = 15
	# calculate PMFs
	masspoints = 0:size
	pmf_a = dcarpbin(masspoints, size=size, shift=shift)
	pmf_b = dao0(masspoints, size=size, steepness=shift, prevalence=0)
	# compare
	round(pmf_b-pmf_a, 3)
})

# test: `rao0` vs `dao0`
set.seed(608)
with(new.env(), {
	# parameters
	sampsize = 4e3
	steepness = -1.22
	size = 10
	prevalence = 0.4
	# generate
	masspoints = 0:size
	x = rao0(sampsize, size=size, steepness=steepness, prevalence=prevalence)
	# calculate PMFs
	pmf_theoretical = setNames(dao0(masspoints, size=size, 
	steepness=steepness, prevalence=prevalence), masspoints)
	pmf_empirical = setNames(sapply(masspoints, function(v) {
		mean(x==v)
	}), masspoints)
	# compare empirical vs theoretical
	plot(masspoints, pmf_empirical, type='h', lwd=2, 
	ylim=c(0, max(pmf_empirical)))
	points(masspoints, pmf_theoretical, pch=1, lwd=2)
	round(pmf_empirical-pmf_theoretical, 3)
})

# moment for AO0
moment_ao0 = function(power, size, steepness, prevalence) {
	masspoints = (0:size)^power
	pmf = dao0(masspoints, size=size, steepness=steepness, 
	prevalence=prevalence)
	expected = sum(masspoints*pmf)
	return(expected)
}

# expectation for AO0
expect_ao0 = function(size, steepness, prevalence) {
	return(moment_ao0(power=1, size=size, steepness=steepness, 
	prevalence=prevalence))
}

# variance for AO0
var_ao0 = function(size, steepness, prevalence) {
	masspoints = 0:size
	pmf = dao0(masspoints, size=size, steepness=steepness, 
	prevalence=prevalence)
	first_moment = sum(masspoints*pmf)
	squared_deviation = (masspoints-first_moment)^2
	return(sum(pmf*squared_deviation))
}

# test: expect_ao0 vs rao0
set.seed(412)
(function(x) {
	print(mean(x))
	print(median(x)) 
	boxplot(x)
})(replicate(400, {
	# parameters
	sampsize = 10e3
	steepness = -1
	size = 10
	prevalence = 0.4
	# generate
	masspoints = 0:size
	x = rao0(sampsize, size=size, steepness=steepness, prevalence=prevalence)
	# calculate mean
	expected_theoretical = expect_ao0(size=size, steepness=steepness, 
	prevalence=prevalence)
	expected_empirical = mean(x)
	# compare empirical vs theoretical
	expected_empirical-expected_theoretical
}))

# test: var_ao0 vs rao0
set.seed(412)
(function(x) {
	print(mean(x))
	print(median(x)) 
	boxplot(x)
})(replicate(400, {
	# parameters
	sampsize = 10e3
	steepness = -2.1
	size = 20
	prevalence = 0.6
	# generate
	masspoints = 0:size
	x = rao0(sampsize, size=size, steepness=steepness, prevalence=prevalence)
	# calculate mean
	var_theoretical = var_ao0(size=size, steepness=steepness, 
	prevalence=prevalence)
	var_empirical = var(x)
	# compare empirical vs theoretical
	var_empirical-var_theoretical
}))

# posterior CNR probability for AO0
calc_postr_cnr_ao0 = function(x, size, steepness, prevalence) {
	likelihood_class0 = dcarpbin(x, size=size, shift=steepness)
	likelihood_class1 = dcarpbin(x, size=size, shift=0)
	postr_numerator = prevalence*likelihood_class1
	postr_denominator = postr_numerator+(1-prevalence)*likelihood_class0
	postr = postr_numerator/postr_denominator
	return(postr)
}

# test: `calc_postr_cnr_ao0` vs `rao0`
set.seed(506)
with(new.env(), {
	# parameters
	sampsize = 10e3
	steepness = -2
	size = 10
	prevalence = 0.3
	# generate
	y = sample(0:1, size=sampsize, prob=c(1-prevalence, prevalence), 
	replace=TRUE)
	x_class0 = rcarpbin(sampsize, size=size, shift=steepness)
	x_class1 = rcarpbin(sampsize, size=size, shift=0)
	x = ifelse(y==1, x_class1, x_class0)
	# calculate probabilities
	freq = table(y, x)
	postr_empirical = freq[2,]/colSums(freq)
	postr_theoretical = calc_postr_cnr_ao0(0:size, size=size, 
	steepness=steepness, prevalence=prevalence)
	# compare empirical vs theoretical
	plot(postr_theoretical, postr_empirical); abline(0:1)
	round(postr_empirical-postr_theoretical, 3)
})

# fit method of moments
fit_methmom = function(success_count, size, steepness_lim) {
	# get empirical moments from data
	emp_mean = mean(success_count)
	emp_var = var(success_count)
	# get class 1 (flat) moments
	flat_mean = size/2
	flat_var = sum((0:size)^2)/(size+1)-flat_mean^2
	# get steepest class 0 mean
	lowest_mean = expect_ao0(size=size, steepness=steepness_lim, 
	prevalence=0)
	# emergency exit for extreme empirical mean
	if(emp_mean>flat_mean) {
		fitted_params = list(prevalence=1, steepness=0)
		return(fitted_params)
	} # too high
	if(emp_mean<lowest_mean) {
		fitted_params = list(prevalence=0, steepness=steepness_lim)
		return(fitted_params)
	} # too low
	# set range for steepness
	steepness_left = steepness_lim
	steepness_right = optimize(f=function(s) {
		implied_mean = expect_ao0(size=size, steepness=s, prevalence=0)
		(implied_mean-emp_mean)^2
	}, lower=steepness_lim, upper=0)$minimum
	# helper functions for fitting
	match_prevalence = function(steepness) {
		matched_class0_mean = expect_ao0(size=size, 
		steepness=steepness, prevalence=0)
		matched_prevalence = (emp_mean-matched_class0_mean)/
		(flat_mean-matched_class0_mean)
		matched_prevalence
	} # for given steepness, find prevalence to match mean
	lossfun = function(steepness) {
		matched_prevalence = match_prevalence(steepness)
		implied_mixture_var = var_ao0(size=size, steepness=steepness, 
		prevalence=matched_prevalence)
		abs(emp_var-implied_mixture_var)
	} # for given steepness, find discrepancy in variance
	# optimization
	fitted_steepness = optimize(f=lossfun, lower=steepness_left, 
	upper=steepness_right)$minimum
	fitted_prevalence = match_prevalence(fitted_steepness)
	# result
	fitted_params = list(prevalence=fitted_prevalence, 
	steepness=fitted_steepness)
	return(fitted_params)
}

# test: fitting
set.seed(516)
with(new.env(), {
	# parameters
	sampsize = 2e3
	size = 200
	steepness = -1.21
	prevalence = 0.4
	# generate
	masspoints = 0:size
	x = rao0(sampsize, size=size, steepness=steepness, 
	prevalence=prevalence)
	# fitting
	fit_methmom(success_count=x, size=size, steepness_lim=-10)
})


# small test on CarpBinTable
with(new.env(), {
	# params
	size = 15
	shift_limit = -5
	num_gridpoints = 30
	# implied params
	masspoints = 0:size
	shifts = seq(from=shift_limit, to=0, length.out=num_gridpoints)
	# construct table
	cbtable = CarpBinTable$new(size=size, shift_limit=shift_limit, 
	num_gridpoints=num_gridpoints)
	# vector of errors
	testpoints = -3.12
	pmf_efficient = as.vector(sapply(testpoints, function(s) {
		cbtable$dcarpbin(shift=s, success_counts=masspoints)
	}))
	pmf_baseline = as.vector(sapply(testpoints, function(s) {
		dcarpbin(masspoints, size=size, shift=s)
	}))
	plot(pmf_efficient, pmf_baseline, main=Sys.time()); abline(0:1)
	round(pmf_efficient-pmf_baseline, 4)
})

# big test on shifts not in table
Sys.time(); with(new.env(), {
	# params
	size = 200
	shift_limit = -10
	num_gridpoints = 200
	# implied params
	masspoints = 0:size
	shifts = seq(from=shift_limit, to=0, length.out=num_gridpoints)
	# construct table
	cbtable = CarpBinTable$new(size=size, shift_limit=shift_limit, 
	num_gridpoints=num_gridpoints)
	# matrix of errors
	testpoints = (shifts[-1]+
	(shifts[-length(shifts)]))/2 # in between gridpoints!
	tab_efficient = as.vector(sapply(testpoints, function(s) {
		cbtable$dcarpbin(shift=s, success_counts=masspoints)
	}))
	tab_baseline = sapply(testpoints, function(s) {
		dcarpbin(masspoints, size=size, shift=s)
	})
	tab_diff = tab_efficient-tab_baseline
	# display result
	if(FALSE) {
		rownames(tab_diff) = masspoints
		colnames(tab_diff) = round(testpoints, 5)
		image(x=rownames(tab_diff), y=colnames(tab_diff), z=tab_diff, 
		main=Sys.time())
	}
	round(quantile(abs(tab_diff)), 3)
}); Sys.time()

# test lookup mean in CarpBinTable
with(new.env(), {
	# params
	size = 200
	shift_limit = -1
	num_gridpoints = 300
	# implied params
	masspoints = 0:size
	shifts = seq(from=shift_limit, to=0, length.out=num_gridpoints)
	# construct table
	cbtable = CarpBinTable$new(size=size, shift_limit=shift_limit, 
	num_gridpoints=num_gridpoints)
	testpoints = (shifts[-1]+
	(shifts[-length(shifts)]))/2 # in between gridpoints!
	# test 1
	baseline = sapply(testpoints, function(s) {
		expect_ao0(size=size, steepness=s, prevalence=0)
	})
	efficient = cbtable$lookup_m(testpoints)
	plot(testpoints, baseline, lwd=3, col='blue', ylab='mean', 
	main=Sys.time())
	lines(testpoints, efficient, lwd=3, col='red')
	err1 = baseline-efficient
	print(quantile(round(abs(err1), 3)))
	# test 2
	testpoints_again = cbtable$reverse_lookup_m(efficient)
	err2 = testpoints-testpoints_again
	print(quantile(round(abs(err2), 3)))
})

# test lookup variance in CarpBinTable
with(new.env(), {
	# params
	size = 200
	shift_limit = -1
	num_gridpoints = 300
	# implied params
	masspoints = 0:size
	shifts = seq(from=shift_limit, to=0, length.out=num_gridpoints)
	# construct table
	cbtable = CarpBinTable$new(size=size, shift_limit=shift_limit, 
	num_gridpoints=num_gridpoints)
	testpoints = (shifts[-1]+
	(shifts[-length(shifts)]))/2 # in between gridpoints!
	# test 1
	baseline = sapply(testpoints, function(s) {
		var_ao0(size=size, steepness=s, prevalence=0)
	})
	efficient = cbtable$lookup_v(testpoints)
	plot(testpoints, sqrt(baseline), lwd=3, col='blue', ylab='stdev', 
	main=Sys.time())
	lines(testpoints, sqrt(efficient), lwd=3, col='red')
	err1 = sqrt(baseline)-sqrt(efficient)
	print(quantile(round(abs(err1), 3)))
	# test 2
	testpoints_again = cbtable$reverse_lookup_v(efficient)
	err2 = testpoints-testpoints_again
	print(quantile(round(abs(err2), 3)))
})

# test: AO1 posterior calculation
with(new.env(), {
	# parameters
	size = 200
	steepness = -1.411
	prevalence = 0.11
	shift_limit = -10
	num_gridpoints = 200
	# create objects
	cbtable = CarpBinTable$new(size=size, shift_limit=shift_limit, 
	num_gridpoints=num_gridpoints)
	mod = AO1Model$new(table=cbtable)
	print(class(mod))
	# calculate probabilities
	masspoints = 0:size
	baseline = calc_postr_cnr_ao0(masspoints, size=size, 
	steepness=steepness, prevalence=prevalence)
	efficient = mod$calc_postr_cnr(success_counts=masspoints, 
	steepness=steepness, features=cbind(rep(1, 
	length.out=length(masspoints))), slopes=rbind(qlogis(prevalence)))
	# compare
	plot(baseline, efficient, main=Sys.time()); abline(0:1)
	err = efficient-baseline
	quantile(round(abs(err), 3))
})

# test: AO0 posterior calculation
with(new.env(), {
	# parameters
	size = 200
	steepness = -1.411
	prevalence = 0.11
	shift_limit = -10
	num_gridpoints = 200
	# create objects
	cbtable = CarpBinTable$new(size=size, shift_limit=shift_limit, 
	num_gridpoints=num_gridpoints)
	mod = AO0Model$new(table=cbtable)
	# calculate probabilities
	masspoints = 0:size
	baseline = calc_postr_cnr_ao0(masspoints, size=size, 
	steepness=steepness, prevalence=prevalence)
	mod$par_set(steepness=steepness, prevalence=prevalence)
	efficient = mod$calc_postr_cnr(success_counts=masspoints)
	# compare
	plot(baseline, efficient, main=Sys.time()); abline(0:1)
	err = efficient-baseline
	quantile(round(abs(err), 3))
})

# test: AO0 PMF calculation
with(new.env(), {
	# parameters
	size = 200
	steepness = -3.11
	prevalence = 0.21
	shift_limit = -10
	num_gridpoints = 200
	# create objects
	cbtable = CarpBinTable$new(size=size, shift_limit=shift_limit, 
	num_gridpoints=num_gridpoints)
	mod = AO0Model$new(tables=cbtable)
	# calculate probabilities
	masspoints = 0:size
	baseline = dao0(masspoints, size=size, steepness=steepness, 
	prevalence=prevalence)
	efficient = mod$dao0(masspoints, steepness=steepness, 
	prevalence=prevalence)
	# compare
	plot(efficient, baseline, main=Sys.time()); abline(0:1)
	err = efficient-baseline
	quantile(round(abs(err), 3))
})

# test: AO0 fitting
set.seed(226)
with(new.env(), {
	# parameters
	sampsize = 100
	size = 200
	steepness = -2.14
	prevalence = 0.4
	shift_limit = -5
	num_gridpoints = 300
	# generate
	y = sample(0:1, size=sampsize, prob=c(1-prevalence, prevalence), 
	replace=TRUE)
	sc_class0 = rcarpbin(sampsize, size=size, shift=steepness)
	sc_class1 = rcarpbin(sampsize, size=size, shift=0)
	sc = ifelse(y==1, sc_class1, sc_class0)
	# create objects
	cbtable = CarpBinTable$new(size=size, shift_limit=shift_limit, 
	num_gridpoints=num_gridpoints)
	mod_ml = AO0Model$new(tables=cbtable)
	mod_mm = AO0Model$new(tables=cbtable)	
	# stuff to display
	display = function(mod) {
		print(mod$loglikelihood())
		print(lapply(mod$par_get()[-1], round, digits=3))
		met = mod$calc_metrics(true_class_labels=y)
		metric_shortnames = setNames(c('confmat', 'acc', 'sens', 
		'spec', 'ppv', 'npv', 'flagrate'), c('confusion', 'accuracy', 
		'sensitivity', 'specificity', 'positive_predictive_value', 
		'negative_predictive_value', 'flag_rate'))
		names(met) = metric_shortnames[names(met)]
		print(round(unlist(met[-1]), 3))
	}
	# fit maximum likelihood
	message('ML')
	mod_ml$fit(sc, init=c(-1, 0.3))
	display(mod_ml)
	# fit method of moments
	message('MM')
	mod_mm$fit_mm(sc)
	display(mod_mm)
})

# test: AO1 fitting
set.seed(221)
with(new.env(), {
	# parameters
	sampsize = 150
	size = 200
	steepness = -1.21
	slopes = c(-0.5*1, c(1, -1, +2, -2))
	features = cbind(1, replicate(length(slopes)-1, {
		runif(sampsize, -1, +1)
	}))
	shift_limit = -5
	num_gridpoints = 300
	# generate
	prevalence = as.vector(plogis(features%*%slopes))
	y = rbinom(length(prevalence), size=1, prob=prevalence)
	sc_class0 = rcarpbin(sampsize, size=size, shift=steepness)
	sc_class1 = rcarpbin(sampsize, size=size, shift=0)
	sc = ifelse(y==1, sc_class1, sc_class0)
	# create objects
	cbtable = CarpBinTable$new(size=size, shift_limit=shift_limit, 
	num_gridpoints=num_gridpoints)
	mod1_ml = AO1Model$new(tables=cbtable)
	mod1_em = AO1Model$new(tables=cbtable)
	mod0 = AO0Model$new(tables=cbtable)
	# stuff to display
	display = function(mod) {
		print(mod$loglikelihood())
		print(lapply(mod$par_get()[-1], round, digits=3))
		met = mod$calc_metrics(true_class_labels=y)
		metric_shortnames = setNames(c('confmat', 'acc', 'sens', 
		'spec', 'ppv', 'npv', 'flagrate'), c('confusion', 'accuracy', 
		'sensitivity', 'specificity', 'positive_predictive_value', 
		'negative_predictive_value', 'flag_rate'))
		names(met) = metric_shortnames[names(met)]
		print(round(unlist(met[-1]), 3))
	}
	# AO0 MM
	message('AO0-MM')
	mod0$fit_mm(sc)
	display(mod0)
	# AO1 EM
	message('EM')
	mmest = mod0$par_get()
	mod1_em$fit_em(sc, features, maxit=10, init=NULL, tol=0.01, 
	verbose=TRUE)
	display(mod1_em)
	# AO1 canned ML
	message('ML')
	mod1_ml$fit(sc, features, init=NULL)
	display(mod1_ml)
})

# test: decision boundary
set.seed(137)
with(new.env(), {
	# parameters
	sampsize = 200
	size = 200
	steepness = -0.81
	features = cbind(1, runif(sampsize, -1, +1))
	slopes = c(-1, 2)
	shift_limit = -5
	num_gridpoints = 300
	# generate
	prevalence = as.vector(plogis(features%*%slopes))
	y = rbinom(length(prevalence), size=1, prob=prevalence)
	print(table(y))
	sc_class0 = rcarpbin(sampsize, size=size, shift=steepness)
	sc_class1 = rcarpbin(sampsize, size=size, shift=0)
	sc = ifelse(y==1, sc_class1, sc_class0)
	# create carpbin tables
	cbtable = CarpBinTable$new(size=size, shift_limit=shift_limit, 
	num_gridpoints=num_gridpoints)
	# stuff to display
	display = function(mod) {
		print(mod$loglikelihood())
		print(lapply(mod$par_get()[-1], round, digits=3))
		met = mod$calc_metrics(true_class_labels=y)
		metric_shortnames = setNames(c('confmat', 'acc', 'sens', 
		'spec', 'ppv', 'npv', 'flagrate'), c('confusion', 'accuracy', 
		'sensitivity', 'specificity', 'positive_predictive_value', 
		'negative_predictive_value', 'flag_rate'))
		names(met) = metric_shortnames[names(met)]
		print(met$confmat)
		print(round(unlist(met[-1]), 3))
	}
	# set true params AO1
	mod1 = AO1Model$new(tables=cbtable)
	mod1$data_set(success_counts=sc, features=features)
	mod1$par_set(steepness=steepness, slopes=slopes)
	# set true params AO0
	mod0 = AO0Model$new(tables=cbtable)
	mod0$data_set(success_counts=sc)
	mod0$par_set(steepness=steepness, prevalence=mean(prevalence))
	# draw decision boundary AO1
	message('AO1')
	display(mod1)
	boundary_coords = mod1$coords_decibo()
	correct1 = y==round(mod1$calc_postr_cnr())
	plot(qlogis(prevalence), sc, pch=as.character(y), 
	col=ifelse(correct1, 'green', 'red'), xlab='covariate', 
	ylab='success count')
	lines(boundary_coords[['lincomb']], boundary_coords[['boundarycount']])
	# draw decision boundary AO0
	message('AO0')
	display(mod0)
	sc_threshold = mod0$threshold_ao0()
	correct0 = y==round(mod0$calc_postr_cnr())
	abline(h=sc_threshold)
})

