# generate from carpal distribution
rcarpal = function(n, shift) {
	z = rnorm(n)
	success_rate = pnorm(z+shift)
	return(success_rate)
}

# carpal density
dcarpal = function(x, shift) {
	in_range = 0<=x&x<=1
	density_in_range = dnorm(qnorm(x)-shift)/dnorm(qnorm(x))
	return(ifelse(in_range, density_in_range, 0))
}

# generate from carpal binomial hierarchy
rcarpbin = function(n, size, shift) {
	success_rate = rcarpal(n, shift=shift)
	success_count = rbinom(n, size=size, prob=success_rate)
	return(success_count)
}

# carpal binomial hierarchy PMF
dcarpbin = function(x, size, shift) {
	# function for single inputs
	f = function(x, size, shift) {
		fun_to_integrate = function(u) {
			#dbinom(x, size=size, prob=u)*dcarpal(u, shift=shift)
			dbinom(x, size=size, prob=pnorm(u))*dnorm(u, mean=shift)
		}
		integrate(fun_to_integrate, lower=-Inf, upper=+Inf)$value
	}
	# apply for multiple inputs
	return(mapply(f, x=x, size=size, shift=shift))
}

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
		stopifnot(length(shift)==1)
		# find index to use
		idx_lower = findInterval(shift, private$shifts, 
		rightmost.closed=TRUE)
		if(idx_lower==0|idx_lower==length(private$shifts)) {
			return(NA)
		} # out of range
		idx_upper = idx_lower+1
		# interpolate
		betweenness = (shift-private$shifts[idx_lower])/
		(private$shifts[idx_upper]-private$shifts[idx_lower])
		pmf_lower = private$table_pmf[idx_lower, ]
		pmf_upper = private$table_pmf[idx_upper, ]
		pmf = pmf_lower+betweenness*(pmf_upper-pmf_lower)
		# look up interpolated values
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
	# set data, features only
	data_set_features = function(features) {
		private$features = features
		return(invisible(NULL))
	},
	# set data, success counts only
	data_set_success_counts = function(success_counts) {
		private$success_counts = success_counts
		return(invisible(NULL))
	},
	# clear data
	data_clear = function() {
		self$data_set(success_counts=NULL, features=NULL)
		return(invisible(NULL))
	},
	# compute prevalences
	calc_prevalence_ao1 = function(features=private$features, 
	slopes=private$slopes) {
		prevalence = as.vector(plogis(features%*%slopes))
		return(prevalence)
	},
	# compute prevalences (again)
	calc_prevalence = function(features=private$features, 
	slopes=private$slopes) {
		return(self$calc_prevalence_ao1(features=features, 
		slopes=slopes))
	},
	# compute class likelihoods
	calc_class_likelihoods = function(success_counts, 
	steepness=private$steepness) {
		pmf_class0 = private$tables$dcarpbin(shifts=steepness, 
		success_counts=success_counts)
		pmf_class1 = private$tables$dcarpbin(shifts=0, 
		success_counts=success_counts)
		df = setNames(data.frame(pmf_class0, pmf_class1), c('0', '1'))
		return(df)
	},
	# compute prior times likelihood
	calc_prior_times_likelihood_ao1 = function(
	success_counts=private$success_counts, features=private$features, 
	steepness=private$steepness, slopes=private$slopes) {
		prev_cnr = self$calc_prevalence(features=features, 
		slopes=slopes)
		prevalences = setNames(list(1-prev_cnr, prev_cnr), c('0', '1'))
		class_likelihoods = self$calc_class_likelihoods(
		success_counts=success_counts, steepness=steepness)
		product = mapply(function(u1, u2) {
			u1*u2
		}, u1=prevalences, u2=class_likelihoods)
		return(product)
	},
	# compute prior times likelihood (again)
	calc_prior_times_likelihood = function(
	success_counts=private$success_counts, features=private$features, 
	steepness=private$steepness, slopes=private$slopes) {
		return(self$calc_prior_times_likelihood_ao1(success_counts, 
		features=features, steepness=steepness, slopes=slopes))
	},
	# compute AO1 PMF
	dao1 = function(success_counts=private$success_counts, 
	features=private$features, steepness=private$steepness, 
	slopes=private$slopes) {
		product = self$calc_prior_times_likelihood_ao1(
		success_counts=success_counts, features=features, 
		steepness=steepness, slopes=slopes)
		prob_mix = rowSums(product)
		return(prob_mix)
	},
	# compute likelihood
	loglikelihood_ao1 = function(steepness=private$steepness, 
	slopes=private$slopes, features=private$features, 
	success_counts=private$success_counts) {
		lik = self$dao1(success_counts=success_counts, features=features,
		steepness=steepness, slopes=slopes)
		return(sum(log(lik)))
	},
	# compute likelihood (again)
	loglikelihood = function(steepness=private$steepness, 
	slopes=private$slopes, features=private$features, 
	success_counts=private$success_counts) {
		return(self$loglikelihood_ao1(steepness=steepness, slopes=slopes, 
		features=features, success_counts=success_counts))
	},
	# compute posterior probability
	calc_postr_ao1 = function(success_counts=private$success_counts, 
	features=private$features, steepness=private$steepness, 
	slopes=private$slopes) {
		product = self$calc_prior_times_likelihood_ao1(
		success_counts=success_counts, features=features, 
		steepness=steepness, slopes=slopes)
		class_postr = t(sapply(1:nrow(product), function(i) {
			nume = product[i, ]
			deno = sum(nume)
			nume/deno
		}))
		return(class_postr)
	},
	# compute posterior probability (again)
	calc_postr = function(success_counts=private$success_counts, 
	features=private$features, steepness=private$steepness, 
	slopes=private$slopes) {
		return(self$calc_postr_ao1(success_counts=success_counts, 
		features=features, steepness=steepness, slopes=slopes))
	},
	# calculate bayes classifier metrics
	calc_metrics_ao1 = function(true_class_labels, 
	success_counts=private$success_counts, steepness=private$steepness, 
	features=private$features, slopes=private$slopes) {
		# get predicted labels
		class1prob = self$calc_postr_ao1(success_counts=success_counts,
		features=features, steepness=steepness, slopes=slopes)[, '1']
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
		return(self$calc_metrics_ao1(true_class_labels, 
		success_counts=success_counts, steepness=steepness, 
		features=features, slopes=slopes))
	},
	# conditional threshold
	calc_boundary = function(linear_comb, steepness=private$steepness) {
		masspoints = 0:private$size
		intercept = cbind(rep(1, times=length(masspoints)))
		boundary = sapply(linear_comb, function(u) {
			# get posterior probability of CNR
			postr_cnr = self$calc_postr_ao1(
			success_counts=masspoints, features=intercept, 
			steepness=steepness, 
			slopes=u)[, '1'] # use AO1 method
			# emergency exit if all predict same class
			round_postr_cnr = round(postr_cnr)
			if(round_postr_cnr[1]==
			round_postr_cnr[length(round_postr_cnr)]) {
				return(NA)
			}
			# interpolate
			idx_lower = max(which(postr_cnr<=0.5))
			idx_upper = min(which(postr_cnr>=0.5))
			idx = c(idx_lower, idx_upper)
			approx(x=postr_cnr[idx], y=masspoints[idx], 
			xout=0.5)$y
		})
		return(boundary)
	},
	
	
	# coordinates for decision boundary
	coords_decibo_ao1 = function(success_counts=private$success_counts,
	steepness=private$steepness, features=private$features, 
	slopes=private$slopes, true_class_labels=NULL) {
		# get all mass points
		masspoints = 0:private$size
		# compute the linear combination values
		lin_comb = as.vector(as.matrix(features)%*%slopes)
		# for each linear combination value, get the boundary count
		boundary_count = sapply(lin_comb, self$calc_boundary)
		# get the predicted class label		
		predicted_class_labels = round(self$calc_postr()[, '1'])
		# make sorted dataframe
		df = data.frame(lincomb=lin_comb, boundarycount=boundary_count, 
		rawcount=success_counts, predlabel=predicted_class_labels)
		if(!is.null(true_class_labels)) {
			df[['truelabel']] = true_class_labels
		}
		df_ordered = df[order(df[['lincomb']]), ]
		return(df_ordered)
	},
	# coordinates for decision boundary (again)
	coords_decibo = function(success_counts=private$success_counts,
	steepness=private$steepness, features=private$features, 
	slopes=private$slopes, true_class_labels=NULL) {
		return(self$coords_decibo_ao1(success_counts=success_counts, 
		steepness=steepness, features=features, slopes=slopes, 
		true_class_labels=true_class_labels))
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
			fitted = list(steepness=0, prevalence=1)
			return(fitted)
		} # too high
		if(emp_mean<lowest_mean) {
			fitted = list(steepness=shift_limit, prevalence=0)
			return(fitted)
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
	}
))

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

	# cast as AO1
	cast_as_ao1 = function(success_counts, prevalence=plogis(private$slopes)) {
		features = cbind(rep(1, times=length(success_counts)))
		slopes = qlogis(prevalence)
		common = list(features=features, slopes=slopes)
		return(common)
	},
	# compute prior times likelihood
	calc_prior_times_likelihood = function(success_counts, 
	steepness=private$steepness, prevalence=plogis(private$slopes)) {
		# cast as AO1
		common = self$cast_as_ao1(success_counts=success_counts, 
		prevalence=prevalence)
		# use AO1 method
		product = self$calc_prior_times_likelihood_ao1(
		success_counts=success_counts, features=common$features, 
		steepness=steepness, slopes=common$slopes)
		return(product)
	},
	# compute AO0 PMF
	dao0 = function(success_counts, steepness=private$steepness, 
	prevalence=plogis(private$slopes)) {
		product = self$calc_prior_times_likelihood(
		success_counts=success_counts, steepness=steepness, 
		prevalence=prevalence)
		prob_mix = rowSums(product)
		return(prob_mix)
	},
	# compute likelihood
	loglikelihood = function(steepness=private$steepness, 
	prevalence=plogis(private$slopes), success_counts=private$success_counts) {
		lik = self$dao0(success_counts=success_counts, steepness=steepness, 
		prevalence=prevalence)
		return(sum(log(lik)))
	},
	# compute posterior probability
	calc_postr = function(success_counts=private$success_counts, 
	steepness=private$steepness, prevalence=plogis(private$slopes)) {
		product = self$calc_prior_times_likelihood(
		success_counts=success_counts, steepness=steepness, 
		prevalence=prevalence)
		class_postr = t(sapply(1:nrow(product), function(i) {
			nume = product[i, ]
			deno = sum(nume)
			nume/deno
		}))
		return(class_postr)
	},
	# calculate bayes classifier metrics
	calc_metrics = function(true_class_labels, 
	success_counts=private$success_counts, steepness=private$steepness, 
	prevalence=plogis(private$slopes)) {
		# cast as AO1
		common = self$cast_as_ao1(success_counts=success_counts, 
		prevalence=prevalence)
		# use AO1 method
		metrics = self$calc_metrics_ao1(true_class_labels, 
		success_counts=success_counts, steepness=steepness, 
		features=common$features, slopes=common$slopes)
		return(metrics)
	},
	# coordinates for decision boundary
	coords_decibo = function(success_counts=private$success_counts,
	steepness=private$steepness, slopes=plogis(private$slopes), 
	true_class_labels=NULL) {
		# cast as AO1
		common = self$cast_as_ao1(success_counts=success_counts, 
		prevalence=prevalence)
		# use AO1 method
		coords = self$coords_decibo_ao1(
		success_counts=success_counts, steepness=steepness, 
		features=common$features, slopes=common$slopes, 
		true_class_labels=true_class_labels)
		return(coords)
	},

	# get AO0 decision boundary
	threshold_ao0 = function(steepness=private$steepness, 
	prevalence=plogis(private$slopes)) {
		return(self$calc_boundary(qlogis(prevalence), 
		steepness=steepness))
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

