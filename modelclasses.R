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
		pmf = private$table[shift_idx, ]
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
	plot(pmf_efficient, pmf_baseline); abline(0:1)
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
		image(x=rownames(tab_diff), y=colnames(tab_diff), z=tab_diff)
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
	plot(testpoints, baseline, lwd=3, col='blue', ylab='mean')
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
	plot(testpoints, sqrt(baseline), lwd=3, col='blue', ylab='stdev')
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
	initialize = function(table) {
		# assert
		stopifnot('CarpBinTable' %in% class(table))
		# use table information
		private$table = table
		private$size = private$tables$par()$size
	},
	# get params
	par = function() {
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
	data = function() {
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
		dcarpbin = private$tables$dcarpbin
		pmf_class0 = dcarpbin(shifts=steepness, 
		success_counts=success_counts)
		pmf_class1 = dcarpbin(shifts=0, 
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
		likelihood_class0 = dcarpbin(shifts=private$steepness, 
		success_counts=success_counts)
		likelihood_class1 = dcarpbin(shifts=0, 
		success_counts=success_counts)
		# numerator and denominator of posterior
		postr_numerator = prevalence*likelihood_class1
		postr_denominator = postr_numerator+(1-prevalence)*likelihood_class0
		postr = postr_numerator/postr_denominator
		return(postr)
	},
	fit_canned = function(success_counts, init) {
		# define objective
		objective = function(u) {
			likelihood = self$dao0(success_counts, steepness=u[1], 
			prevalence=u[2])
			-1*sum(log(likelihood))
		}
		# optimize
		shift_limit = private$table$par()$shift_limit
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
		private$steepness = estimate[1]
		private$prevalence = estimate[2]
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
	mod = AO0Model$new(table=cbtable, prevalence=prevalence, steepness=steepness)
	# calculate probabilities
	masspoints = 0:size
	postr_a = calc_postr_cnr_ao0(masspoints, size=size, 
	steepness=steepness, prevalence=prevalence) # non-R6 implementation
	postr_b = mod$calc_postr_cnr(masspoints) # R6 implementation
	# compare non-R6 implementation vs R6
	plot(postr_a, postr_b); abline(0:1)
	round(postr_a-postr_b, 3)
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
	mod = AO0Model$new(table=cbtable, prevalence=prevalence, 
	steepness=steepness)
	# calculate probabilities
	masspoints = 0:size
	pmf_a = dao0(masspoints, size=size, steepness=steepness, 
	prevalence=prevalence) # non-R6 implementation
	pmf_b = mod$dao0(masspoints, steepness=steepness, 
	prevalence=prevalence) # R6 implementation
	# compare non-R6 implementation vs R6
	plot(pmf_a, pmf_b); abline(0:1)
	round(pmf_a-pmf_b, 3)
})

# test: fitting
set.seed(516)
with(new.env(), {
	# parameters
	sampsize = 2e3
	size = 200
	steepness = -2.14
	prevalence = 0.4
	shift_limit = -5
	num_gridpoints = 300
	# generate
	masspoints = 0:size
	x = rao0(sampsize, size=size, steepness=steepness, prevalence=prevalence)
	# create objects
	cbtable = CarpBinTable$new(size=size, shift_limit=shift_limit, 
	num_gridpoints=num_gridpoints)
	mod = AO0Model$new(table=cbtable, prevalence=NA, steepness=NA)
	print(c(steepness, prevalence))
	mod$fit_canned(x, init=c(-1, 0.5))
	unlist(mod$par())
})
