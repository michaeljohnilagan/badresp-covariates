# load analysis results
load('./simulation.RData')

# enumerateo n x contam pairs
eg = expand.grid(n=sim_factors$n, contam=sim_factors$contam)

# range of accuracies
print(range(as.matrix(sim_tab[, c('met_sc.acc', 'met_ao0.acc', 
'met_ao1.acc')])))

# rename n and contam
sim_tab_panelplots = sim_tab
sim_tab_panelplots$n = paste('n=', sim_tab$n, sep='')
sim_tab_panelplots$contam = paste('contam=', sim_tab$contam, sep='')

# do panel plots: accuracy
pdf('figs-simulation-acc.pdf')
classifier2col = setNames(c('red3', 'green3', 'blue3'), c('sc', 'ao0', 'ao1'))
key_acc = list(columns=3, lines=list(lty=rep(1:2, times=3), 
col=rep(classifier2col[c('sc', 'ao0', 'ao1')], each=2)), 
text=list(label=c('all items, SC', 'even items only, SC', 'all items, AO0', 
'even items only, AO0', 'all items, AO1', 'even AO1')))
lattice::xyplot(met_sc.acc~factor(xsep, levels=sim_factors$xsep) | 
factor(n)*factor(contam), data=sim_tab_panelplots, panel=function(x,y,...) {
	# get panel data
	i = as.numeric(lattice::panel.number()[1])
	curr_ncontam = subset(sim_tab, n==eg$n[i]&contam==eg$contam[i])
	# plot for zwhich=='all'
	coords_all = curr_ncontam[curr_ncontam$zwhich=='all', ]
	x_all = factor(coords_all$xsep, levels=sim_factors$xsep)
	y_all_sc = coords_all$met_sc.acc # SC
	lattice::panel.lines(x_all, y_all_sc, lty=1, col=classifier2col['sc'], 
	type='b', pch=19)
	y_all_ao0 = coords_all$met_ao0.acc # AO0
	lattice::panel.lines(x_all, y_all_ao0, lty=1, col=classifier2col['ao0'], 
	type='b', pch=19)
	y_all_ao1 = coords_all$met_ao1.acc # AO1
	lattice::panel.lines(x_all, y_all_ao1, lty=1, col=classifier2col['ao1'], 
	type='b', pch=19)
	# plot for zwhich=='even'
	coords_even = curr_ncontam[curr_ncontam$zwhich=='even', ]
	x_even = factor(coords_even$xsep, levels=sim_factors$xsep)
	y_even_sc = coords_even$met_sc.acc # SC 
	lattice::panel.lines(x_even, y_even_sc, lty=2, col=classifier2col['sc'], 
	type='b', pch=19)
	y_even_ao0 = coords_even$met_ao0.acc # AO0
	lattice::panel.lines(x_even, y_even_ao0, lty=2, col=classifier2col['ao0'], 
	type='b', pch=19)
	y_even_ao1 = coords_even$met_ao1.acc # AO1
	lattice::panel.lines(x_even, y_even_ao1, lty=2, col=classifier2col['ao1'], 
	type='b', pch=19)
}, ylim=c(0.6, 1), xlab='separation in covariates', ylab='accuracy', 
key=key_acc)
dev.off()

# do panel plots: LL improved
pdf('figs-simulation-llimproved.pdf')
key_ll = list(lines=list(lty=1:2), text=list(label=c('all', 'even')))
lattice::xyplot(ll_improved~factor(xsep, levels=sim_factors$xsep) | 
factor(n)*factor(contam), data=sim_tab_panelplots, panel=function(x,y,...) {
	# get panel data
	i = as.numeric(lattice::panel.number()[1])
	curr_ncontam = subset(sim_tab, n==eg$n[i]&contam==eg$contam[i])
	# plot for zwhich=='all'
	coords_all = curr_ncontam[curr_ncontam$zwhich=='all', ]
	x_all = factor(coords_all$xsep, levels=sim_factors$xsep)
	y_all = coords_all$ll_improved
	lattice::panel.lines(x_all, y_all, lty=1, col=1, type='b', pch=19)
	# plot for zwhich=='even'
	coords_even = curr_ncontam[curr_ncontam$zwhich=='even', ]
	x_even = factor(coords_even$xsep, levels=sim_factors$xsep)
	y_even = coords_even$ll_improved 
	lattice::panel.lines(x_even, y_even, lty=2, col=1, type='b', pch=19)
}, xlab='separation in covariates', ylab='rate of log likelihood improved',
key=key_ll)
dev.off()

# end session
Sys.time()
devtools::session_info()
