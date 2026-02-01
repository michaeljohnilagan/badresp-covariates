# load analysis results
load('./analysis.RData')

# create figs directory if not exists
(function(path) {
	if(!dir.exists(path)) {
		dir.create(path)
	}
})('./figs')

# function: plot the unsupervised result, for a pointscale
visualize_ps_unsuperv = function(true_class_labels, model_ao0, model_ao1, 
nomsens, main=NA) {
	# hard coded graphical params
	lwd = 3
	cex.axis = 2
	cex.lab = 1.5
	cex = 1.5
	# colors for correct vs incorrect prediction
	correct_col = 'gray40'
	incorrect_col = 'red'
	# get coordinates
	coords = model_ao1$coords_decibo(true_class_labels=true_class_labels)
	# plot AO1 decision boundary
	with(coords, {
		plot(lincomb, rawcount, type='n', 
		xlab='linear combination of covariates', ylab='success count', 
		main=main, 
		cex=cex, cex.axis=cex.axis, cex.lab=cex.lab)
		lines(lincomb, boundarycount, lwd=lwd)
	})
	# plot AO0 decision boundary
	thresh_ao0 = model_ao0$threshold_ao0()
	abline(h=thresh_ao0, lty=1, lwd=lwd)
	# plot sensitivity calibrated decision boundary
	size = model_ao1$par_get()$size
	thresh_sc = (1-nomsens)*(size+1)-1
	abline(h=thresh_sc, lty=2, lwd=lwd)
	# plot points
	with(subset(coords, truelabel==predlabel), {
		points(lincomb, rawcount, pch=as.character(truelabel), 
		col=correct_col, cex=cex)
	}) # correct predictions
	with(subset(coords, truelabel!=predlabel), {
		points(lincomb, rawcount, pch=as.character(truelabel), 
		col=incorrect_col, cex=cex)
	}) # incorrect predictions
	return(NULL)
}

# produce PDF for unsupervised decision boundaries
pdf('./figs/figs-analysis-decibo.pdf')
invisible(sapply(results_unsuperv, function(o) {
	# hard coded graphical params
	lwd = 3
	cex = 1.5
	bg = 'gray90'
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
	# legend
	if(FALSE) {
		legend('topleft', title='boundaries', lty=3:1,
		legend=c('AO1', 'AO0', '95% SC'), bg=bg, cex=cex, lwd=lwd)
		legend('left', title='result', pch=19, lty=0, 
		col=c('gray40', 'red'), legend=c('correct', 'incorrect'), 
		bg=bg, cex=cex, lwd=lwd)
	}
}))
dev.off()

# produce PDF for unsupervised decision boundaries
pdf('./figs/figs-analysis-boxplot.pdf')
invisible(sapply(results_unsuperv, function(o) {
	# hard coded graphical params
	cex.axis = 2
	cex.lab = 1.5
	cex = 1.5
	pointscale = o$plotting$pointscale
	# unpack
	y = o$plotting$y
	pval = o$plotting$fit$pval
	# plot
	boxplot(pval~y, ylab='p value', xlab='true class label', 
	main=paste(pointscale, ' point scale', sep=''),
	cex=cex, cex.axis=cex.axis, cex.lab=cex.lab)
}))
dev.off()

# produce PDF for unsupervised decision boundaries
pdf('./figs/figs-analysis-ecdf.pdf')
invisible(sapply(results_unsuperv, function(o) {
	# hard coded graphical params
	lwd = 2
	cex.axis = 2
	cex.lab = 1.5
	cex = 1.35
	pointscale = o$plotting$pointscale
	bg = 'gray90'
	# unpack
	y = o$plotting$y
	pval = o$plotting$fit$pval
	# plot
	cdf0 = ecdf(pval[y==0])
	cdf1 = ecdf(pval[y==1])
	plot(cdf0, col='blue', ylim=c(0, 1), xlim=c(0, 1),
	xlab='p value', ylab='quantile rank', 
	main=paste(pointscale, ' point scale', sep=''), 
	cex=cex, cex.axis=cex.axis, cex.lab=cex.lab, lwd=lwd)
	lines(cdf1, col='red', lwd=lwd, cex=cex)
	abline(0:1, col='black', lty=1, lwd=lwd)
	legend('bottomright', legend=c('not careless', 'careless'),
	pch=19, lty=1, col=c('blue', 'red'), bg=bg, cex=cex)
}))
dev.off()

# function: visualize unsupervised vs supervised
visualize_compare_by_outcome_measure = function(tables, outcome_measure, 
name2col, ylab=NA) {
	# assert
	stopifnot(length(outcome_measure)==1)
	stopifnot(outcome_measure %in% c('acc', 'sens', 'spec', 'ppv', 'npv', 
	'flagrate'))
	# hard coded graphical params
	lwd = 3
	cex.axis = 2
	cex.lab = 1.5
	cex = 1.5
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
	xlab=xlab, ylab=ylab,
	cex.axis=cex.axis, cex.lab=cex.lab, cex=cex)
	axis(1, at=pointscales, labels=pointscales, cex.axis=cex.axis)
	# determine colors
	name2col = name2col
	names(name2col) = paste(names(name2col), '.', outcome_measure, sep='')
	# draw lines
	column_names = colnames(relevant)
	sapply(column_names, function(nm) {
		lines(pointscales, relevant[, nm], col=name2col[nm], type='b', 
		lwd=lwd, pch=19)
	})
	# output the numbers
	return(relevant)
}

# produce PDF for unsupervised vs supervised
pdf('./figs/figs-analysis-compare.pdf')
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
	legend=name2legend[classifier_names], cex=1.1, bg='gray90')
	# plot sensitivity
	vis_sens = visualize_compare_by_outcome_measure(tables=relevant_metrics, 
	outcome_measure='sens', ylab='sensitivity', name2col=name2col)
	legend('right', pch=19, col=name2col[classifier_names], 
	legend=name2legend[classifier_names], cex=1.1, bg='gray90')
	# plot specificity
	vis_spec = visualize_compare_by_outcome_measure(tables=relevant_metrics, 
	outcome_measure='spec', ylab='specificity', name2col=name2col)
	legend('right', pch=19, col=name2col[classifier_names], 
	legend=name2legend[classifier_names], cex=1.1, bg='gray90')
})
dev.off()

# end session
Sys.time()
devtools::session_info()
