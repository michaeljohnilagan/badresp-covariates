source('funs.R')

# efficient table lookup class
CarpBinTable = R6::R6Class('CarpBinTable', 
private=list(
	size = NULL,
	shifts = NULL,
	table_pmf = NULL,
	table_mv = NULL
), public=list(
	# initialize
	initialize = function(size, shift_limit, num_gridpoints) {
		# assert
		stopifnot(shift_limit<0) # negative shifts only
		# set size
		private$size = size
		# set shifts
		private$shifts = seq(from=shift_limit, to=0, 
		length.out=num_gridpoints)
		# set PMF table
		masspoints = 0:size
		private$table_pmf = sapply(masspoints, function(k) {
			sapply(private$shifts, function(s) {
				dcarpbin(k, size=size, shift=s)
			}) # rows are shifts
		}) # columns are mass points
		# set mean variance table
		private$table_mv = data.frame(shift=private$shifts)
		private$table_mv[['mean']] = private$table_pmf%*%masspoints
		private$table_mv[['variance']] = sapply(
		1:nrow(private$table_mv), function(i) {
			masspoints_sqdev = (masspoints-
			private$table_mv[['mean']][i])^2
			sum(private$table_pmf[i, ]*masspoints_sqdev)
		})
	},
	# get params
	par = function() {
		list_of_params = list(size=private$size, 
		shift_limit=private$shifts[1], 
		num_gridpoints=length(private$shifts))
		return(list_of_params)
	},
	# compute dcarpbin for one shift value
	dcarpbin_sameshift = function(shift, success_counts) {
		# assert
		stopifnot(shift<=0)
		stopifnot(length(shift)==1)
		# look up
		shift_idx = min(which(shift<=private$shifts))
		pmf = private$table_pmf[shift_idx, ]
		masspoints = 0:private$size
		lookedup = pmf[match(success_counts, masspoints)]
		return(lookedup)
	},
	# compute dcarpbin for multiple shift values
	dcarpbin = function(shifts, success_counts) {
		if(length(shifts)==1) {
			# case: only one shift value
			lookedup = self$dcarpbin_sameshift(shift=shifts, 
			success_counts=success_counts)
		} else {
			# case: multiple shift values
			framed = data.frame(shifts=shifts, 
			success_counts=success_counts)
			lookedup = sapply(1:nrow(framed), function(i) {
				shift = framed[['shifts']][i]
				success_count = framed[['success_counts']][i]
				self$dcarpbin_sameshift(shift=shift, 
				success_counts=success_count)
			})
		}
		return(lookedup)
	},
	# given shift, lookup mean
	lookup_m = function(shift) {
		approx(x=private$table_mv[['shift']], 
		y=private$table_mv[['mean']], xout=shift)$y
	},
	# given mean, lookup shift
	reverse_lookup_m = function(mean) {
		approx(x=private$table_mv[['mean']], 
		y=private$table_mv[['shift']], xout=mean)$y
	},
	# given shift, lookup variance
	lookup_v = function(shift) {
		approx(x=private$table_mv[['shift']], 
		y=private$table_mv[['variance']], xout=shift)$y
	},
	# given variance, lookup shift
	reverse_lookup_v = function(var) {
		approx(x=private$table_mv[['variance']], 
		y=private$table_mv[['shift']], xout=var)$y
	}
))

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

# model class AO0
AO0Model = R6::R6Class('AO0Model', 
private=list(
	size = NULL,
	steepness = NULL,
	prevalence = NULL,
	success_counts = NULL,
	tables = NULL
), public=list(
	# initialize object
	initialize = function(tables) {
		# assert
		stopifnot('CarpBinTable' %in% class(tables))
		# use table information
		private$tables = tables
		private$size = private$tables$par()$size
	},
	# get params
	par_get = function() {
		list_of_params = list(size=private$size, 
		steepness=private$steepness, prevalence=private$prevalence)
		return(list_of_params)
	},
	# set params
	par_set = function(steepness=NULL, prevalence=NULL) {
		private$steepness = steepness
		private$prevalence = prevalence
		return(invisible(NULL))
	},
	# clear params
	par_clear = function() {
		self$par_set(steepness=NULL, prevalence=NULL)
		return(invisible(NULL))
	},
	# get data
	data_get = function() {
		return(private$success_counts)
	},
	# set data
	data_set = function(success_counts) {
		private$success_counts = success_counts
		return(invisible(NULL))
	},
	# clear data
	data_clear = function() {
		self$data_set(success_counts=NULL)
		return(invisible(NULL))
	},
	# compute PMF
	dao0 = function(success_counts, steepness, prevalence) {
		# compute PMF by class
		pmf_class0 = private$tables$dcarpbin(shifts=steepness, 
		success_counts=success_counts)
		pmf_class1 = private$tables$dcarpbin(shifts=0, 
		success_counts=success_counts)
		# mix the two classes
		pmf_mix = prevalence*pmf_class1+(1-prevalence)*pmf_class0
		return(pmf_mix)
	},
	# compute posterior probability of CNR
	calc_postr_cnr = function(success_counts=private$success_counts, 
	steepness=private$steepness, prevalence=private$prevalence) {
		# compute likelihood by class
		dcarpbin = private$tables$dcarpbin
		likelihood_class0 = dcarpbin(shifts=steepness, 
		success_counts=success_counts)
		likelihood_class1 = dcarpbin(shifts=0, 
		success_counts=success_counts)
		# numerator and denominator of posterior
		postr_numerator = prevalence*likelihood_class1
		postr_denominator = postr_numerator+(1-prevalence)*
		likelihood_class0
		postr = postr_numerator/postr_denominator
		return(postr)
	},
	# fitting via maximum likelihood
	fit_ml = function(success_counts, init) {
		# define objective
		objective = function(u) {
			lik = self$dao0(success_counts=success_counts, 
			steepness=u[1], prevalence=u[2])
			-1*sum(log(lik))
		}
		# optimize
		shift_limit = private$tables$par()$shift_limit
		if(is.null(init)) {
			init = c(shift_limit/2, 0.5)
		}
		if(TRUE) {
			estimate = optim(init, fn=objective)$par # using nelder-mead
		} else {
			estimate = optim(init, fn=objective, 
			method='L-BFGS-B', lower=c(shift_limit, 0), 
			upper=c(0, 1))$par # using L-BFGS-B
		}
		# assign result to object
		self$par_set(steepness=estimate[1], prevalence=estimate[2])
		self$data_set(success_counts=success_counts)
		return(invisible(NULL))
	},
	# find params implied by matching to mean for fixed steepness
	match_implied_params = function(steepness, mix_mean) {
		flat_mean = private$size/2
		class0_mean = private$tables$lookup_m(steepness)
		prevalence = (mix_mean-class0_mean)/(flat_mean-class0_mean)
		class0_var = private$tables$lookup_v(steepness)
		matched = list(prevalence=prevalence, class0_mean=class0_mean, 
		class0_var=class0_var)
		return(matched)
	},
	# fitting via method of moments
	fit_mm = function(success_counts) {
		# get empirical moments from data
		emp_mean = mean(success_counts)
		emp_var = var(success_counts)
		# get class 1 (flat) moments
		flat_mean = private$size/2
		flat_var = sum((0:private$size)^2)/
		(private$size+1)-flat_mean^2
		# get steepest class 0 mean
		shift_limit = private$tables$par()$shift_limit
		lowest_mean = private$tables$lookup_m(shift_limit)
		# emergency exit for extreme empirical mean
		if(emp_mean>flat_mean) {
			self$par_set(steepness=0, prevalence=1)
			self$data_set(success_counts=success_counts)
			return(invisible(NULL))
		} # too high
		if(emp_mean<lowest_mean) {
			self$par_set(steepness=shift_limit, prevalence=0)
			self$data_set(success_counts=success_counts)
			return(invisible(NULL))
		} # too low
		# set range for steepness
		steepness_left = shift_limit
		steepness_right = private$tables$reverse_lookup_m(emp_mean)
		# loss function for each steepness
		lossfun = function(steepness) {
			# find matching prevalence and class 0 params
			matched = self$match_implied_params(steepness, mix_mean=emp_mean)
			# compute implied mixture variance
			mse_class0_part = matched[['class0_var']]+
			(matched[['class0_mean']]-emp_mean)^2
			mse_class1_part = flat_var+(flat_mean-emp_mean)^2
			implied_mix_var = (1-matched[['prevalence']])*
			mse_class0_part+matched[['prevalence']]*mse_class1_part
			# compare
			abs(implied_mix_var-emp_var)
		}
		# optimization
		fitted_steepness = optimize(f=lossfun, lower=steepness_left, 
		upper=steepness_right)$minimum
		fitted_prevalence = self$match_implied_params(fitted_steepness, 
		mix_mean=emp_mean)[['prevalence']]
		# assign result to object
		self$par_set(steepness=fitted_steepness, prevalence=fitted_prevalence)
		self$data_set(success_counts=success_counts)
		return(invisible(NULL))
	}
))

# test: posterior calculation
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
	efficient = mod$calc_postr_cnr(success_counts=masspoints, 
	steepness=steepness, prevalence=prevalence)
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

# test: fitting
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
	x_class0 = rcarpbin(sampsize, size=size, shift=steepness)
	x_class1 = rcarpbin(sampsize, size=size, shift=0)
	x = ifelse(y==1, x_class1, x_class0)
	# create objects
	cbtable = CarpBinTable$new(size=size, shift_limit=shift_limit, 
	num_gridpoints=num_gridpoints)
	mod = AO0Model$new(tables=cbtable)
	print(c(steepness, prevalence))
	# fit maximum likelihood
	mod$fit_ml(x, init=c(-1, 0.3))
	print(unlist(mod$par_get())[-1])
	postr_ml = mod$calc_postr_cnr()
	predbin_ml = round(postr_ml)
	acc_ml = mean(predbin_ml==y)
	print(acc_ml)
	# fit method of moments
	mod$fit_mm(x)
	print(unlist(mod$par_get())[-1])
	postr_mm = mod$calc_postr_cnr()
	predbin_mm = round(postr_mm)
	acc_mm = mean(predbin_mm==y)
	print(acc_mm)
})

# model class AO1
AO1Model = R6::R6Class('AO1Model', 
private=list(
	size = NULL,
	steepness = NULL,
	slopes = NULL,
	features = NULL,
	success_counts = NULL,
	tables = NULL
), public=list(
	# initialize object
	initialize = function(tables) {
		# assert
		stopifnot('CarpBinTable' %in% class(tables))
		# use table information
		private$tables = tables
		private$size = private$tables$par()$size
	},
	# get params
	par_get = function() {
		list_of_params = list(size=private$size, 
		steepness=private$steepness, slopes=private$slopes)
		return(list_of_params)
	},
	# set params
	par_set = function(steepness=NULL, slopes=NULL) {
		private$steepness = steepness
		private$slopes = slopes
		return(invisible(NULL))
	},
	# clear params
	par_clear = function() {
		self$par_set(steepness=NULL, slopes=NULL)
		return(invisible(NULL))
	},
	# get data
	data_get = function() {
		list_of_data = list(success_counts=private$success_counts,
		features=private$features)
		return(list_of_data)
	},
	# set data
	data_set = function(success_counts, features) {
		private$success_counts = success_counts
		private$features = features
		return(invisible(NULL))
	},
	# clear data
	data_clear = function() {
		self$data_set(success_counts=NULL, features=NULL)
		return(invisible(NULL))
	},
	# compute prevalences
	calc_prev = function(features, slopes) {
		prevalences = plogis(features%*%slopes)
		return(prevalences)
	},
	# compute AO0 PMF
	dao0 = function(success_counts, steepness, prevalences) {
		# compute PMF by class
		pmf_class0 = private$tables$dcarpbin(shifts=steepness, 
		success_counts=success_counts)
		pmf_class1 = private$tables$dcarpbin(shifts=0, 
		success_counts=success_counts)
		# mix the two classes
		pmf_mix = prevalences*pmf_class1+(1-prevalences)*pmf_class0
		return(pmf_mix)
	},
	# compute AO1 PMF
	dao1 = function(success_counts, steepness, features, slopes) {
		prevalences = self$calc_prev(features=features, slopes=slopes)
		pmf_mix = self$dao0(success_counts=success_counts, 
		steepness=steepness, prevalences=prevalences)
		return(pmf_mix)
	},
	# compute posterior probability of CNR
	calc_postr_cnr = function(success_counts=private$success_counts, 
	steepness=private$steepness, features=private$features, 
	slopes=private$slopes) {
		# compute class 1 (CNR) prior
		prevalence = self$calc_prev(features=features, slopes=slopes)
		# compute likelihood by class
		likelihood_class0 = private$tables$dcarpbin(shifts=steepness, 
		success_counts=success_counts)
		likelihood_class1 = private$tables$dcarpbin(shifts=0, 
		success_counts=success_counts)
		# numerator and denominator of posterior
		postr_numerator = prevalence*likelihood_class1
		postr_denominator = postr_numerator+(1-prevalence)*
		likelihood_class0
		postr = postr_numerator/postr_denominator
		return(postr)
	},
	# fitting via maximum likelihood
	fit_ml = function(success_counts, features, init=NULL) {
		# define objective as negative log likelihood
		negloglik = function(u) {
			lik = self$dao1(success_counts=success_counts, 
			steepness=u[1], features=features, slopes=u[-1])
			-1*sum(log(lik))
		}
		# optimize
		shift_limit = private$tables$par()$shift_limit
		if(is.null(init)) {
			init = c(shift_limit/2, rep(0, times=ncol(features)))
		}
		estimate = optim(init, fn=negloglik)$par # using nelder-mead
		# force valid steepness
		if(estimate[1]>0) {
			estimate[1] = 0
		}
		# assign result to object
		self$par_set(steepness=estimate[1], slopes=estimate[-1])
		self$data_set(success_counts=success_counts, features=features)
		return(invisible(NULL))
	}
))

# test: posterior calculation
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
	steepness=steepness, features=cbind(rep(1, length.out=length(masspoints))), 
	slopes=rbind(qlogis(prevalence)))
	# compare
	plot(baseline, efficient, main=Sys.time()); abline(0:1)
	err = efficient-baseline
	quantile(round(abs(err), 3))
})

# test: fitting
set.seed(226)
with(new.env(), {
	# parameters
	sampsize = 300
	size = 200
	steepness = -2.14
	prevalence = 0.4
	shift_limit = -5
	num_gridpoints = 300
	# generate
	y = sample(0:1, size=sampsize, prob=c(1-prevalence, prevalence), 
	replace=TRUE)
	x_class0 = rcarpbin(sampsize, size=size, shift=steepness)
	x_class1 = rcarpbin(sampsize, size=size, shift=0)
	x = ifelse(y==1, x_class1, x_class0)
	# create objects
	cbtable = CarpBinTable$new(size=size, shift_limit=shift_limit, 
	num_gridpoints=num_gridpoints)
	mod = AO1Model$new(tables=cbtable)
	print(c(steepness, qlogis(prevalence)))
	# fit maximum likelihood
	mod$fit_ml(x, cbind(rep(1, times=length(x))), init=NULL)
	print(unlist(mod$par_get())[-1])
	postr_ml = mod$calc_postr_cnr()
	predbin_ml = round(postr_ml)
	acc_ml = mean(predbin_ml==y)
	print(acc_ml)
})
