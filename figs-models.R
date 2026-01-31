# import my stuff
Sys.time(); ao = new.env()
with(ao, {
	source('./utils/aomodels.R', local=TRUE)
	source('./utils/classification.R', local=TRUE)
}); Sys.time()

# plot carpaldens
pdf('./figs/figs-models-carpal.pdf')
with(new.env(), {
	# common params
	lwd = 3
	cex.axis = 2
	cex.lab = 1.5
	cex = 1.5
	# draw curves
	curve(ao$dcarpal(x, shift=-2), from=0, to=1, 
	ylab='density', xlab='success rate',
	lty=3, lwd=lwd, cex.axis=cex.axis, cex.lab=cex.lab)
	curve(ao$dcarpal(x, shift=-1), from=0, to=1, add=TRUE, lty=2, lwd=lwd)
	curve(ao$dcarpal(x, shift=0), from=0, to=1, add=TRUE, lty=1, lwd=lwd)
	# put legend
	legend('topright', lty=3:1, cex=cex, lwd=lwd,
	legend=c(expression(delta==-2), expression(delta==-1), 
	expression(delta==0)))
})
dev.off()

# function for barplots
ao0_barplot = function(steepness, prevalence, size, max_height=NULL, 
header=TRUE, shift_limit=-10, num_gridpoints=300) {
	# graphical params
	if(is.null(max_height)) {
		ylim = NULL
	} else {
		ylim = c(0, max_height)
	}
	cex.names = 1.0
	cex.axis = 1.35
	cex.lab = 1.5
	# initialize model object
	tables_carpbin = ao$CarpBinTable$new(size=size, shift_limit=shift_limit, 
	num_gridpoints=num_gridpoints)
	mod0 = ao$AO0Model$new(tables=tables_carpbin)
	# set AO0 params
	mod0$par_set(steepness=steepness, prevalence=prevalence)
	# calculate matrix for bar plot
	masspoints = 0:size
	pmf_by_class = mod0$calc_prior_times_likelihood(masspoints)
	rownames(pmf_by_class) = masspoints
	# plotting
	if(header) {
		main = bquote(paste(lambda==.(prevalence), ', ', 
		delta==.(steepness)))
	} else {
		main = NULL
	}
	barplot(t(pmf_by_class[, c('1', '0')]), col=c('red', 'blue'), ylim=ylim,
	xlab='success count', ylab='probability', main=main,
	cex.names=cex.names, cex.axis=cex.axis, cex.lab=cex.lab)
	return(pmf_by_class)
}

# AO0 PMF-by-class barplots
pdf("./figs/figs-models-ao0.pdf")
ao0_barplot(prevalence=0.25, steepness=-1, size=15, max_height=0.5, 
header=TRUE)
legend("topright", legend=c('not careless', 'careless'), col=c('blue', 'red'), 
pch=16, cex=1.35, bg="gray90")
ao0_barplot(prevalence=0.25, steepness=-2, size=15, max_height=0.5,
header=TRUE)
ao0_barplot(prevalence=0.75, steepness=-1, size=15, max_height=0.5,
header=TRUE)
ao0_barplot(prevalence=0.75, steepness=-2, size=15, max_height=0.5,
header=TRUE)
dev.off()

# generate data to use in demo for decision boundary
set.seed(2246)
env_decibo = new.env()
Sys.time(); with(env_decibo, {
	# set nominal sensitivity
	nomsens = 0.95
	# initialize model object
	num_perms = 200
	shift_limit = -10
	num_gridpoints = 300
	tables_carpbin = ao$CarpBinTable$new(size=num_perms, shift_limit=shift_limit, 
	num_gridpoints=num_gridpoints)
	mod1 = ao$AO1Model$new(tables=tables_carpbin)
	mod0 = ao$AO0Model$new(tables=tables_carpbin)
	# generate and set features
	n = 300
	covariates_mu = -1
	covariates_sig = 1
	features = cbind(1, rnorm(n, covariates_mu, covariates_sig))
	mod1$data_set_features(features=features)
	# set slopes and steepness
	slopes = c(-1, 2)
	steepness = -2
	mod1$par_set(steepness=steepness, slopes=slopes)
	prevalence = mod1$calc_prevalence(features=features, slopes=slopes)
	mod0$par_set(steepness=steepness, prevalence=mean(prevalence))
	# generate labels
	true_label = rbinom(n, size=1, prob=prevalence)
	# generate success counts
	success_count = ifelse(true_label==1, ao$rcarpbin(n, num_perms, 0), 
	ao$rcarpbin(n, num_perms, steepness))
	pval = (1+success_count)/(1+num_perms)
	mod1$data_set_success_counts(success_counts=success_count)
	mod0$data_set(success_counts=success_count)
	# make coordinates
	coords = mod1$coords_decibo(true_class_labels=true_label)
	coords[['covariate']] = (coords$lincomb-slopes[1])/slopes[2]
}); Sys.time()

# plot demo of decision boundary
pdf('./figs/figs-models-ao1.pdf')
with(env_decibo, {
	# graphical params
	lwd = 3
	cex.axis = 2
	cex.lab = 1.5
	cex = 1.5
	xlim = NULL
	# canvas
	with(coords, {
		plot(covariate, rawcount, type='n', xlim=xlim,
		xlab='covariate', ylab='success count',
		cex.axis=cex.axis, cex.lab=cex.lab, cex=cex)
	})
	# draw boundaries
	with(coords, {
		lines(covariate, boundarycount, lwd=lwd) # AO1 boundary
		thresh_ao0 = mod0$threshold_ao0()
		abline(h=thresh_ao0, lty=1, lwd=lwd) # AO1 boundary
		thresh_sc = (1-nomsens)*(num_perms+1)-1
		abline(h=thresh_sc, lty=2, lwd=lwd) # SC boundary
	})
	# draw points
	correct = with(coords, predlabel==truelabel)
	correct_col = 'gray40'
	incorrect_col = 'red'
	with(subset(coords, correct), {
		points(covariate, rawcount, pch=as.character(truelabel), 
	 	col=correct_col, cex=cex)
	}) # correct points
	with(subset(coords, !correct), {
		points(covariate, rawcount, pch=as.character(truelabel), 
		col=incorrect_col, cex=cex)
	}) # incorrect points
	# legend
	legend('topleft', title='boundaries', lty=3:1,
	legend=c('AO1', 'AO0', '95% SC'), bg='gray', cex=cex, lwd=lwd,)
	legend('left', title='result', pch=19, lty=0, 
	col=c(correct_col, incorrect_col),
	legend=c('correct', 'incorrect'), bg='gray', cex=cex, lwd=lwd)
})
dev.off()

# difference in adjacent PMF
pdf("./figs/figs-models-adjpmfdiff.pdf")
curve(dbinom(5, size=10, prob=x)-dbinom(6, size=10, prob=x), from=0, to=1,
xlab="success rate", 
ylab="probability of 5 minus probability of 6 successes", 
lwd=2, cex=1.5, cex.axis=1.5, cex.lab=1.5);
abline(h=0, lty=3)
dev.off()
