source('funs.R')

# interpolator class
CarpBinInterpolator = R6::R6Class('CarpBinInterpolator', 
private=list(
	size = NULL,
	shifts = NULL,
	table = NULL
), public=list(
	initialize = function(size, shift_limit, num_gridpoints) {
		private$size = size
		stopifnot(shift_limit<0) # negative shifts only
		private$shifts = seq(from=shift_limit, to=0, 
		length.out=num_gridpoints)
		private$table = sapply(0:size, function(k) {
			sapply(private$shifts, function(s) {
				dcarpbin(k, size=size, shift=s)
			}) # rows are shifts
		}) # columns are mass points
	},
	dcarpbin = function(shift, success_counts) {
		# assert
		stopifnot(shift<=0)
		stopifnot(length(shift)==1)
		# take the relevant columns of the lookup table
		masspoints = 0:private$size
		relevant_masspoints = intersect(success_counts, masspoints)
		relevant_tab = private$table[, 1+relevant_masspoints]
		# interpolate probabilities using the relevant subtable
		relevant_probs = sapply(1:ncol(relevant_tab), 
		function(j) {
			approx(x=private$shifts, y=relevant_tab[, j], xout=shift)$y
		})
		# lookup probabilities for the counts of interest
		lookedup = relevant_probs[match(success_counts, 
		relevant_masspoints)]
		return(lookedup)
	},
	par = function() {
		return(list(size=private$size, shift_limit=private$shifts[1], 
		num_gridpoints=length(private$shifts)))
	}
))

# test interpolator on public and private
with(new.env(), {
	# params
	size = 15
	shift_limit = -5
	num_gridpoints = 50
	# construct interpolator
	interpolator = CarpBinInterpolator$new(size=size, shift_limit=shift_limit, 
	num_gridpoints=num_gridpoints)
	interpolator$par()
})

# small test on interpolator
with(new.env(), {
	# params
	size = 15
	shift_limit = -5
	num_gridpoints = 50
	# implied params
	masspoints = 0:size
	shifts = seq(from=shift_limit, to=0, length.out=num_gridpoints)
	# construct interpolator
	interpolator = CarpBinInterpolator$new(size=size, shift_limit=shift_limit, 
	num_gridpoints=num_gridpoints)
	# vector of errors
	testpoints = -4.9
	pmf_efficient = as.vector(sapply(testpoints, function(s) {
		interpolator$dcarpbin(shift=s, success_counts=masspoints)
	}))
	pmf_baseline = as.vector(sapply(testpoints, function(s) {
		dcarpbin(masspoints, size=size, shift=s)
	}))
	plot(pmf_efficient, pmf_baseline); abline(0:1)
	round(pmf_efficient-pmf_baseline, 4)
})

# big test on interpolator
Sys.time(); with(new.env(), {
	# params
	size = 200
	shift_limit = -10
	num_gridpoints = 200
	# implied params
	masspoints = 0:size
	shifts = seq(from=shift_limit, to=0, length.out=num_gridpoints)
	# construct interpolator
	interpolator = CarpBinInterpolator$new(size=size, shift_limit=shift_limit, 
	num_gridpoints=num_gridpoints)
	# matrix of errors
	testpoints = (shifts[-1]+
	(shifts[-length(shifts)]))/2 # in between gridpoints!
	tab_efficient = as.vector(sapply(testpoints, function(s) {
		interpolator$dcarpbin(shift=s, success_counts=masspoints)
	}))
	tab_baseline = sapply(testpoints, function(s) {
		dcarpbin(masspoints, size=size, shift=s)
	})
	tab_diff = tab_efficient-tab_baseline
	# display result
	rownames(tab_diff) = masspoints
	colnames(tab_diff) = round(testpoints, 5)
	if(FALSE) {
		image(x=rownames(tab_diff), y=colnames(tab_diff), z=tab_diff)
	}
	quantile(abs(tab_diff))
}); Sys.time()

# model class AO0
AccOptim = R6::R6Class('AccOptim', 
private=list(
	steepness = NA,
	prevalence = NA,
	interpolator = NULL
), public=list(
	initialize = function(interpolator=NULL, size=NULL, shift_limit=NULL, 
	num_gridpoints=NULL) {
		# use provided interpolator, otherwise use other args
		private$interpolator = if(is.null(interpolator)) {
			CarpBinInterpolator$new(size=size, shift_limit=shift_limit, 
			num_gridpoints=num_gridpoints)
		} else {
			interpolator
		}
	},
	par = function() {
		params = list(size=unname(private$interpolator$parameters()['size']), 
		steepness=private$steepness, prevalence=private$prevalence)
		return(unlist(params))
	}
))
