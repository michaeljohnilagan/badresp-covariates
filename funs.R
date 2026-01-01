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

# convert p value to success count
pval2count = function(pval, size, tolerance=1e-5) {
	success_count = pval*(size+1)-1
	rounding_error = round(success_count)-success_count
	if(max(abs(rounding_error))>tolerance) {
		warning(c('rounding error up to ', 
		rounding_error))
	}
	return(round(success_count))
}

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

