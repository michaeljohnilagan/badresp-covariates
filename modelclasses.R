source('funs.R', local=TRUE)

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
		prevalences = as.vector(plogis(features%*%slopes))
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
	dao1 = function(success_counts, features, steepness, slopes) {
		prevalences = self$calc_prev(features=features, slopes=slopes)
		pmf_mix = self$dao0(success_counts=success_counts, 
		steepness=steepness, prevalences=prevalences)
		return(pmf_mix)
	},
	# compute likelihood
	loglikelihood = function(steepness=private$steepness, 
	slopes=private$slopes, features=private$features, 
	success_counts=private$success_counts) {
		lik = self$dao1(success_counts=success_counts, features=features,
		steepness=steepness, slopes=slopes)
		return(sum(log(lik)))
	},
	# compute posterior probability of CNR
	calc_postr_cnr_ao1 = function(success_counts=private$success_counts, 
	features=private$features, steepness=private$steepness, 
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
	# compute posterior probability of CNR (again)
	calc_postr_cnr = function(success_counts=private$success_counts, 
	features=private$features, steepness=private$steepness, 
	slopes=private$slopes) {
		self$calc_postr_cnr_ao1(success_counts=success_counts, 
		features=features, steepness=steepness, slopes=slopes)
	},
	# calculate bayes classifier metrics
	calc_metrics_ao1 = function(true_class_labels, 
	success_counts=private$success_counts, steepness=private$steepness, 
	features=private$features, slopes=private$slopes) {
		# get predicted labels
		class1prob = self$calc_postr_cnr_ao1(success_counts=success_counts,
		features=features, steepness=steepness, slopes=slopes)
		predicted_class_labels = round(class1prob)
		# metrics
		confmat = table(pred=predicted_class_labels, 
		true=true_class_labels)
		acc = mean(true_class_labels==predicted_class_labels)
		sens = mean(predicted_class_labels[true_class_labels==1]==1)
		spec = mean(predicted_class_labels[true_class_labels==0]==0)
		ppv = mean(true_class_labels[predicted_class_labels==1]==1)
		npv = mean(true_class_labels[predicted_class_labels==0]==0)
		flagrate = mean(predicted_class_labels==1)
		# put together
		metrics = list(confusion=confmat, accuracy=acc, 
		sensitivity=sens, specificity=spec, 
		positive_predictive_value=ppv, negative_predictive_value=npv,
		flag_rate=flagrate)
		return(metrics)
	},
	# calculate bayes classifier metrics (again)
	calc_metrics = function(true_class_labels, 
	success_counts=private$success_counts, steepness=private$steepness, 
	features=private$features, slopes=private$slopes) {
		self$calc_metrics_ao1(true_class_labels=true_class_labels, 
		success_counts=success_counts, steepness=steepness, 
		features=features, slopes=slopes)
	},
	# coordinates for decision boundary
	coords_decibo = function(success_counts=private$success_counts,
	steepness=private$steepness, features=private$features, 
	slopes=private$slopes) {
		# get all mass points
		masspoints = 0:private$size
		# compute the linear combination values
		lin_comb = as.vector(as.matrix(features)%*%slopes)
		# for each linear combination value, get the boundary count
		boundary_count = sapply(lin_comb, function(u) {
			postr = self$calc_postr_cnr(
			success_counts=masspoints, 
			features=cbind(rep(1, times=length(masspoints))),
			slopes=u)
			approx(x=postr, y=masspoints, xout=0.5)$y
		})
		# make sorted dataframe
		df = data.frame(lincomb=lin_comb, boundarycount=boundary_count, 
		rawcount=success_counts)
		df_ordered = df[order(df[['lincomb']]), ]
		return(df_ordered)
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
	match_moments = function(success_counts) {
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
			matched = self$match_implied_params(steepness, 
			mix_mean=emp_mean)
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
		fitted_prevalence = self$match_implied_params(
		fitted_steepness, mix_mean=emp_mean)[['prevalence']]
		# return result
		fitted = list(steepness=fitted_steepness, 
		prevalence=fitted_prevalence)
		return(fitted)
	},
	# fitting via maximum likelihood
	fit = function(success_counts, features, init=NULL) {
		# define objective as negative log likelihood
		negloglik = function(u) {
			loglik = self$loglikelihood(
			success_counts=success_counts, 
			features=as.matrix(features), steepness=u[1], 
			slopes=u[-1])
			-1*loglik
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
		self$data_set(success_counts=success_counts, 
		features=features)
		return(invisible(NULL))
	},
	# fitting via EM
	fit_em = function(success_counts, features, maxit, init=NULL, 
	tol=0.01, verbose=FALSE) {
		# initial values
		if(is.null(init)) {
			init_steepness = self$match_moments(
			success_counts=success_counts)$steepness
			init_slopes = suppressWarnings(cor(success_counts, features))
			init_slopes = ifelse(is.na(init_slopes), 0, init_slopes)
			init = c(init_steepness, init_slopes)
		}
		# projection matrix for least squares
		features = as.matrix(features)
		projection_matrix = solve(t(features)%*%features)%*%t(features)
		# main loop
		old_pars = init
		old_ll = self$loglikelihood(steepness=old_pars[1], 
		slopes=old_pars[-1], features=features, 
		success_counts=success_counts)
		if(verbose) {
			message(c('init | LL: ', round(old_ll, 3)))
		}
		for(iter in 1:maxit) {
			# from old params, calculate posterior
			postr = self$calc_postr_cnr(
			success_counts=success_counts, features=features, 
			steepness=old_pars[1], slopes=old_pars[-1])
			# from posterior, estimate new steepness
			weights_for_moment = (1-postr)/sum(1-postr)
			if(all(weights_for_moment==0)) {
				break
			}
			class0_mean = sum(weights_for_moment*success_counts)
			if(is.na(class0_mean)) {
				break
			}
			new_steepness = private$tables$reverse_lookup_m(
			class0_mean)
			# from posterior, estimate new slopes
			logits = qlogis(postr)
			new_slopes = as.vector(projection_matrix%*%logits)
			# save new candidate
			new_pars = c(new_steepness, new_slopes)
			new_ll = self$loglikelihood(steepness=new_pars[1], 
			slopes=new_pars[-1], features=features, 
			success_counts=success_counts)
			if(verbose) {
				message(c('iter ', iter, ' | LL: ', 
				round(new_ll, 3), ' | mean postr: ', 
				round(mean(postr), 3)))
			}
			# convergence criterion
			converged = sum(abs(new_pars-old_pars))<
			(tol*length(new_pars))
			if(converged) {
				if(verbose) {
					message('converged!')
				}
				break
			}
			# emergency exit if likelihood worsens
			if(new_ll<old_ll|is.na(new_ll)) {
				if(verbose) {
					message('likelihood worsened')
				}
				new_pars = old_pars
				break
			}
			# prepare next iteration
			old_pars = new_pars
			old_ll = new_ll
		}
		# assign result to object
		self$par_set(steepness=new_pars[1], slopes=new_pars[-1])
		self$data_set(success_counts=success_counts, features=features)
		return(invisible(NULL))
	}
))

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

# model class AO0
AO0Model = R6::R6Class('AO0Model', 
inherit=AO1Model,
public=list(
	# get params
	par_get = function() {
		list_of_params = list(size=private$size, 
		steepness=private$steepness, prevalence=plogis(private$slopes))
		return(list_of_params)
	},
	# get params in terms of AO0
	par_get_ao1 = function() {
		list_of_params = list(size=private$size, 
		steepness=private$steepness, slopes=private$slopes)
		return(list_of_params)
	},
	# set params
	par_set = function(steepness=NULL, prevalence=NULL) {
		stopifnot(length(prevalence)==1)
		stopifnot(prevalence>=0&prevalence<=1)
		private$steepness = steepness
		private$slopes = qlogis(prevalence)
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
	# compute likelihood
	loglikelihood = function(steepness=private$steepness, 
	prevalence=plogis(private$slopes), 
	success_counts=private$success_counts) {
		lik = self$dao0(success_counts=success_counts, steepness=steepness, 
		prevalences=prevalence)
		return(sum(log(lik)))
	},
	# compute posterior probability of CNR
	calc_postr_cnr = function(success_counts=private$success_counts, 
	steepness=private$steepness, prevalence=plogis(private$slopes)) {
		features = cbind(rep(1, times=length(success_counts)))
		slopes = qlogis(prevalence)
		self$calc_postr_cnr_ao1(success_counts=success_counts, features=features, 
		steepness=steepness, slopes=slopes)
	},
	# get AO0 decision boundary
	threshold_ao0 = function(steepness=private$steepness, 
	prevalence=plogis(private$slopes)) {
		masspoints = 0:private$size
		postr = self$calc_postr_cnr(success_counts=masspoints, 
		steepness=steepness, prevalence=prevalence)
		threshold = approx(x=postr, y=masspoints, xout=0.5)$y
		return(threshold)
	},
	# calculate bayes classifier metrics
	calc_metrics = function(true_class_labels, 
	success_counts=private$success_counts, steepness=private$steepness, 
	prevalence=plogis(private$slopes)) {
		features = cbind(rep(1, times=length(success_counts)))
		slopes = qlogis(prevalence)
		self$calc_metrics_ao1(true_class_labels=true_class_labels, 
		success_counts=success_counts, features=features, 
		steepness=steepness, slopes=slopes)
	},
	# fitting via method of moments
	fit_mm = function(success_counts) {
		# optimize
		fitted_mm = self$match_moments(success_counts)
		# assign result to object
		self$par_set(steepness=fitted_mm$steepness, 
		prevalence=fitted_mm$prevalence)
		self$data_set(success_counts=success_counts)
		return(invisible(NULL))
	},
	# fitting via maximum likelihood
	fit = function(success_counts, init=NULL) {
		# define objective as negative log likelihood
		negloglik = function(u) {
			lik = self$dao0(success_counts=success_counts, 
			steepness=u[1], prevalence=u[-1])
			-1*sum(log(lik))
		}
		# optimize
		shift_limit = private$tables$par()$shift_limit
		if(is.null(init)) {
			init = c(shift_limit/2, 0.5)
		}
		estimate = optim(init, fn=negloglik)$par # using nelder-mead
		# force valid steepness
		if(estimate[1]>0) {
			estimate[1] = 0
		}
		# assign result to object
		self$par_set(steepness=estimate[1], prevalence=estimate[-1])
		self$data_set(success_counts=success_counts)
		return(invisible(NULL))
	}
))

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

