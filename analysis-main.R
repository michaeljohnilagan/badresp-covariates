# import my stuff
Sys.time(); ao = new.env()
with(ao, {
	source('./utils/aomodels.R', local=TRUE)
	source('./utils/classification.R', local=TRUE)
}); Sys.time()

# load preprocessed data
all_ps = readRDS('./data/analysis-preprocess.RDS')

# initialize objects
settings = list()

# settings: L1P1
settings$l1p1 = with(new.env(), {
	feat_funs = c('mahal', 'ptcor') # nonresponsivity indices
	feat_idvals = c(0, +1) # ideal point
	numperms = 200 # number of synthetic rows per observed
	nomsens = 0.95 # nominal sensitivity
	list(feat_funs=feat_funs, feat_idvals=feat_idvals, numperms=numperms, 
	nomsens=nomsens)
})

# settings: AO models
Sys.time(); settings$ao = with(new.env(), {
	numperms = settings$l1p1$numperms
	shift_limit = -10 # most negative shift possible to fit
	num_gridpoints = 300 # number of grid points to compute
	cbt = ao$CarpBinTable$new(size=numperms, 
	shift_limit=shift_limit, num_gridpoints=num_gridpoints) # table set
	em_maxit = 50 # maximum iterations for EM algorithm
	list(cbt=cbt, em_maxit=em_maxit)
}); Sys.time()

# settings: supervised predictions
settings$superv = with(new.env(), {
	# logistic regression model
	logistic_regression = list(fit=function(df_train) {
		glm(y~., data=df_train, family=binomial(link='logit'))
	}, predict = function(fit, df_test) {
		predict(fit, newdata=df_test, type='response')
	})
	# random forest model
	random_forest = list(fit=function(df_train) {
		randomForest::randomForest(factor(y)~., data=df_train)
	}, predict = function(fit, df_test) {
		predict(fit, newdata=df_test, type='prob')[, '1']
	})
	# number of replicates in cross validation
	num_replicates = 20
	# put together
	list(models=list(lr=logistic_regression, rf=random_forest), 
	num_replicates=num_replicates)
})

# function: impute covariates
impute_covariates = function(test, train) {
	# assert
	stopifnot(ncol(test)==ncol(train))
	# replacement
	test_new = setNames(lapply(1:ncol(test), function(j) {
		current_col_train = train[,j]
		replacement_value = mean(current_col_train[is.finite(current_col_train)], 
		na.rm=TRUE)
		current_col_test = test[,j]
		needs_replacement = is.na(current_col_test)|!is.finite(current_col_test)
		replace(current_col_test, needs_replacement, replacement_value)
	}), colnames(test))
	return(as.data.frame(test_new))
}

# function: run unsupervised analysis for a pointscale
run_ps_unsuperv = function(ps_list, pointscale) {
	# use settings for L1P1
	z = ps_list$z
	pointscales = rep(pointscale, times=ncol(z))
	numperms = settings$l1p1$numperms
	feat_funs = settings$l1p1$feat_funs
	feat_idvals = settings$l1p1$feat_idvals
	nomsens = settings$l1p1$nomsens
	# fit models
	x = ps_list$x
	fit = analysis_unsuperv(z=z, x=x, pointscales=pointscales, 
	numperms=numperms, feat_funs=feat_funs, feat_idvals=feat_idvals)
	# extract predictions
	pred_class_lab_sc = ifelse(fit$pval>=(1-nomsens), 1, 0)
	pred_class_lab_ao0 = round(fit$ao0$calc_postr_cnr())
	pred_class_lab_ao1 = round(fit$ao1$calc_postr_cnr())
	list_preds = list(sc=pred_class_lab_sc, 
	ao0=pred_class_lab_ao0, ao1=pred_class_lab_ao1)
	# calculate metrics
	y = ps_list$y
	metrics = function(predicted_class_label) {
		ao$metrics(true_class_label=y, predicted_class_label=predicted_class_label)
	}
	met = unlist(lapply(list_preds, metrics))
	# stuff for plotting
	plotting = list(y=y, fit=fit, nomsens=nomsens, pointscale=pointscale)
	# put together
	together = list(met=met, plotting=plotting)
	return(together)
}

# function: do unsupervised analysis
analysis_unsuperv = function(z, x, pointscales, numperms, feat_funs, 
feat_idvals, use_em=FALSE) {
	# assert
	stopifnot(nrow(z)==nrow(x))
	stopifnot(length(pointscales)==ncol(z))
	# calculate p values
	p = detranli::cnrdetect(z, pointscales=pointscales, 
	numperms=numperms, feat_funs=feat_funs, feat_idvals=feat_idvals)
	p = ifelse(is.na(p), 1, p)
	succ = ao$pval2count(p, size=numperms)
	# do imputation on features, then add intercept
	x = impute_covariates(test=x, train=x)
	x = cbind(1, as.matrix(x))
	# fit A0 clasiffier
	mod0 = ao$AO0Model$new(tables=settings$ao$cbt)
	mod0$fit_mm(success_counts=succ)
	# fit AO1 classifier
	init = c(mod0$par_get()$steepness, qlogis(mod0$par_get()$prevalence), 
	rep(0, times=ncol(x)-1))
	mod1 = ao$AO1Model$new(tables=settings$ao$cbt)
	if(use_em) {
		trymod1 = try(mod1$fit_em(success_counts=succ, features=x, 
		init=init,maxit=settings$ao$em_maxit, 
		verbose=verbose), silent=TRUE) # EM algorithm
	} else {
		trymod1 = try(mod1$fit(success_counts=succ, features=x, 
		init=init), silent=TRUE) # canned ML
	}
	# special handling if AO1 fails
	optim_exception = class(trymod1)=='try-error'
	if(optim_exception) {
		message('optim exception')
		mod1 = mod0$clone()
	}
	fit = list(pval=p, succ=succ, ao1=mod1, ao0=mod0)
	return(fit)
}

# function: do unsupervised analysis, one replicate of cross validation
analysis_superv = function(id, y, x, num_folds, type) {
	# assert
	stopifnot(length(id)==length(y))
	stopifnot(length(id)==nrow(x))
	# assign folds
	new_order = y+runif(length(y), min=0, max=1)
	id_ordered = id[order(new_order)]
	id2fold = setNames(1+(1:length(id_ordered)) %% num_folds, id_ordered)
	fold = id2fold[id]
	# impute covariates
	x = impute_covariates(test=x, train=x)
	# fit models
	dat_full = as.data.frame(cbind(y, x))
	mods = lapply(1:num_folds, function(k) {
		dat_train = subset(dat_full, fold!=k)
		mod = settings$superv$models[[type]]$fit(df_train=dat_train)
		mod
	})
	# extract predictions
	yhat = rep(NA, times=length(y))
	for(k in 1:num_folds) {
		yhat[fold==k] = settings$superv$models[[type]]$predict(fit=mods[[k]], 
		df_test=subset(as.data.frame(x), fold==k))
	}
	# calculate metrics
	met = ao$metrics(true_class_label=y, predicted_class_label=round(yhat))
	names(met) = paste(type, '_', sprintf('%02d', num_folds), 
	'.', names(met), sep='')
	return(met)
}

# function: do unsupervised analysis
run_ps_superv = function(ps_list, num_folds, type, num_replicates) {
	# unpack list
	id = ps_list$id
	x = ps_list$x
	y = ps_list$y
	# repeatedly analyze
	many_repls = replicate(num_replicates, {
		analysis_superv(id=id, y=y, x=x, type=type, num_folds=num_folds)
	})
	return(rowMeans(many_repls))
}

# work unsupervised
set.seed(1144)
Sys.time(); results_unsuperv = mapply(run_ps_unsuperv, ps_list=all_ps, 
pointscale=c(3, 4, 5, 6, 7, 10, 11), SIMPLIFY=FALSE); Sys.time()
metrics_unsuperv = t(sapply(results_unsuperv, function(o) {o$met}))
round(100*metrics_unsuperv, 1)

# work supervised, logistic regression
set.seed(2246)
Sys.time(); metrics_superv_lr_10 = t(mapply(run_ps_superv, ps_list=all_ps, 
num_folds=10, type='lr', 
num_replicates=settings$superv$num_replicates)); Sys.time()
round(100*metrics_superv_lr_10, 1)

# work supervised, random forest
set.seed(2246)
Sys.time(); metrics_superv_rf_10 = t(mapply(run_ps_superv, ps_list=all_ps, 
num_folds=10, type='rf', 
num_replicates=settings$superv$num_replicates)); Sys.time()
round(100*metrics_superv_rf_10, 1)

# end session
Sys.time()
save.image("./analysis.RData")
devtools::session_info()
