# load analysis results
load('./analysis.RData')

# function: plot the unsupervised result, for a pointscale
visualize_ps_unsuperv = function(true_class_labels, model_ao0, model_ao1, 
nomsens, main=NA) {
	# get coordinates
	coords = model_ao1$coords_decibo(true_class_labels=true_class_labels)
	# plot AO1 decision boundary
	with(coords, {
		plot(lincomb, rawcount, type='n', 
		xlab='linear combination of covariates', ylab='success count', main=main)
		lines(lincomb, boundarycount)
	})
	# plot AO0 decision boundary
	thresh_ao0 = model_ao0$threshold_ao0()
	abline(h=thresh_ao0, lty=1)
	# plot sensitivity calibrated decision boundary
	size = model_ao1$par_get()$size
	thresh_sc = (1-nomsens)*(size+1)-1
	abline(h=thresh_sc, lty=2)
	# plot points
	with(subset(coords, truelabel==predlabel), {
		points(lincomb, rawcount, pch=as.character(truelabel))
	}) # correct predictions
	incorrect_col = 'red'
	with(subset(coords, truelabel!=predlabel), {
		points(lincomb, rawcount, pch=as.character(truelabel), col=incorrect_col)
	}) # incorrect predictions
	return(NULL)
}

# produce PDF for unsupervised decision boundaries
pdf('figs-analysis-decibo.pdf')
invisible(sapply(results_unsuperv, function(o) {
	# unpack
	y = o$plotting$y
	model_ao0 = o$plotting$fit$ao0
	model_ao1 = o$plotting$fit$ao1
	nomsens = o$plotting$nomsens
	pointscale = o$plotting$pointscale
	# plot
	visualize_ps_unsuperv(true_class_labels=y, model_ao0=model_ao0, 
	model_ao1=model_ao1, nomsens=nomsens, 
	main=paste(pointscale, ' point scale', sep=''))
}))
dev.off()

# function: visualize unsupervised vs supervised
visualize_compare_by_outcome_measure = function(tables, outcome_measure, 
name2col, ylab=NA) {
	# assert
	stopifnot(length(outcome_measure)==1)
	stopifnot(outcome_measure %in% c('acc', 'sens', 'spec', 'ppv', 'npv', 
	'flagrate'))
	# hard coded point scales
	pointscales = c(3:7, 10:11)
	xlab = 'number of response categories'
	# data to plot
	all_tables = do.call(cbind, tables)
	pattern = paste('.', outcome_measure, sep='')
	relevant = all_tables[, grepl(pattern, colnames(all_tables))]
	# canvas
	plot(NA, NA, type='n', xaxt='n',
	xlim=range(pointscales), ylim=range(relevant),
	xlab=xlab, ylab=ylab)
	axis(1, at=pointscales, labels=pointscales)
	# determine colors
	name2col = name2col
	names(name2col) = paste(names(name2col), '.', outcome_measure, sep='')
	# draw lines
	column_names = colnames(relevant)
	sapply(column_names, function(nm) {
		lines(pointscales, relevant[, nm], col=name2col[nm], type='b', lwd=2, 
		pch=19)
	})
	# output the numbers
	return(relevant)
}

# produce PDF for unsupervised vs supervised
pdf('figs-analysis-compare.pdf')
with(new.env(), {
	# set colors and legend label
	classifier_names = c('sc', 'ao0', 'ao1', 'lr_10', 'rf_10')
	name2col = setNames(c('red3', 'green3', 'blue3', 'orange', 'brown'), 
	classifier_names)
	name2legend = setNames(c('95% sensitivity calibrated', 'AO0', 'AO1', 
	'logistic regression', 'random forest'), classifier_names)
	# gather relevant metrics
	relevant_metrics = list(metrics_unsuperv, metrics_superv_lr_10, 
	metrics_superv_rf_10)
	# plot accuracy
	vis_acc = visualize_compare_by_outcome_measure(tables=relevant_metrics, 
	outcome_measure='acc', ylab='accuracy', name2col=name2col)
	legend('right', pch=19, col=name2col[classifier_names], 
	legend=name2legend[classifier_names])
	# plot sensitivity
	vis_sens = visualize_compare_by_outcome_measure(tables=relevant_metrics, 
	outcome_measure='sens', ylab='sensitivity', name2col=name2col)
	legend('right', pch=19, col=name2col[classifier_names], 
	legend=name2legend[classifier_names])
	# plot specificity
	vis_spec = visualize_compare_by_outcome_measure(tables=relevant_metrics, 
	outcome_measure='spec', ylab='specificity', name2col=name2col)
	legend('right', pch=19, col=name2col[classifier_names], 
	legend=name2legend[classifier_names])
})
dev.off()

# end session
Sys.time()
devtools::session_info()
