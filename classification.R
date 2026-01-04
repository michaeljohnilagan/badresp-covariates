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

# get classification metrics
metrics = function(true_class_labels, predicted_class_labels, 
short_names=FALSE) {
	# compute accuracy, etc
	acc = mean(true_class_labels==predicted_class_labels)
	sens = mean(predicted_class_labels[true_class_labels==1]==1)
	spec = mean(predicted_class_labels[true_class_labels==0]==0)
	ppv = mean(true_class_labels[predicted_class_labels==1]==1)
	npv = mean(true_class_labels[predicted_class_labels==0]==0)
	flagrate = mean(predicted_class_labels==1)
	# put together
	metrics = unlist(list(acc=acc, sens=sens, spec=spec, ppv=ppv, 
	npv=npv, flagrate=flagrate))
	return(metrics)
}

