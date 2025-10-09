source('funs.R')

# AO0 lookup class
ao0lookup = R6::R6Class('ao0lookup', public=list(
	size = NULL,
	shifts = NULL,
	table = NULL,
	initialize = function(size, shift_limit, num_gridpoints) {
		self$size = size
		stopifnot(shift_limit<0)
		self$shifts = seq(from=shift_limit, to=0, length.out=num_gridpoints)
		self$table = sapply(0:size, function(k) {
			sapply(self$shifts, function(s) {
				dcarpbin(k, size=size, shift=s)
			})
		})
	},
	lookup = function(shift, success_counts) {
		# assert
		stopifnot(length(shift)==1)
		stopifnot(shift<=0)
		# work
		masspoints = 0:self$size
		relevant_masspoints = intersect(success_counts, masspoints)
		relevant_pmf = sapply(1+relevant_masspoints, function(j) {
			approx(x=self$shifts, y=self$table[, j], xout=shift)$y
		})
		lookedup = relevant_pmf[match(success_counts, relevant_masspoints)]
		return(lookedup)
	}
))

# small test on lookup
with(new.env(), {
	# params
	size = 15
	shift_limit = -5
	num_gridpoints = 50
	# implied params
	masspoints = 0:size
	shifts = seq(from=shift_limit, to=0, length.out=num_gridpoints)
	# construct interpolator
	interpolator = ao0lookup$new(size=size, shift_limit=shift_limit, 
	num_gridpoints=num_gridpoints)
	# vector of errors
	testpoints = -4.9
	tab_efficient = as.vector(sapply(testpoints, function(s) {
		interpolator$lookup(shift=s, success_counts=masspoints)
	}))
	tab_baseline = as.vector(sapply(testpoints, function(s) {
		dcarpbin(masspoints, size=size, shift=s)
	}))
	plot(tab_efficient, tab_baseline); abline(0:1)
	round(tab_efficient-tab_baseline, 4)
})

# big test on lookup
Sys.time(); with(new.env(), {
	# params
	size = 200
	shift_limit = -10
	num_gridpoints = 200
	# implied params
	masspoints = 0:size
	shifts = seq(from=shift_limit, to=0, length.out=num_gridpoints)
	# construct interpolator
	interpolator = ao0lookup$new(size=size, shift_limit=shift_limit, 
	num_gridpoints=num_gridpoints)
	# matrix of errors
	testpoints = (shifts[-1]+
	(shifts[-length(shifts)]))/2 # in between gridpoints!
	tab_efficient = as.vector(sapply(testpoints, function(s) {
		interpolator$lookup(shift=s, success_counts=masspoints)
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
