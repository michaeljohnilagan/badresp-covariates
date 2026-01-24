# list files
list.files()

# load raw dataframes
env_raw = new.env()
with(env_raw, {
	# likert data
	dat_z = readRDS('./data/experiments1&4-results-to-analysis.RDS')
	# covariates data
	dat_x = read.csv('./data/auxiliary-data-wide.csv', sep=',')
})
sapply(env_raw, dim)
colnames(env_raw$dat_z)
str(env_raw$dat_x)

# verify that rows have unique ID
sapply(env_raw, function(df) {
	all(table(df[['token']])==1)
})

# verify that one is a subset of the other
setdiff(env_raw$dat_x[['token']], env_raw$dat_z[['token']])

# verify that common columns match
with(new.env(), {
	dat_x = env_raw$dat_x
	dat_z = env_raw$dat_z
	common_columns = setdiff(intersect(colnames(dat_x), colnames(dat_z)), 
	'token')
	setNames(lapply(common_columns, function(nm) {
		df = merge(dat_x[c('token', nm)], dat_z[c('token', nm)], 
		by='token')
		table(df[, 2], df[, 3])
	}), common_columns)
})

# group related columns
colnames_list = list(id='token', pointscale='responseScale', y='stimuli', 
z=colnames(env_raw$dat_z)[grepl('_code', colnames(env_raw$dat_z))],
x=c('sumScreenTimeLog', 'engagementScore', 'attentivenessScore', 
'log_dX_relAll', 'sqrt_flipsXAll', 'log_dY_relAll', 'sqrt_flipsYAll', 
'log_vXAll', 'log_vYAll', 'log_aXAll', 'log_aYAll'))
colnames_list

# join data frames
df_unsplit = with(new.env(), {
	colnames_dat_z = unlist(colnames_list[c('id', 'pointscale', 'y', 'z')])
	colnames_dat_x = unlist(colnames_list[c('id', 'x')])
	dat_z = env_raw$dat_z[, colnames_dat_z]
	dat_x = env_raw$dat_x[, colnames_dat_x]
	key = colnames_list$id
	merge(dat_z, dat_x, by=key)
}) # x_from_z not included here
str(df_unsplit)

# check missingness of likert
sort(sapply(unlist(colnames_list$z), function(nm) {
	mean(is.na(df_unsplit[[nm]]))
}), decreasing=TRUE)

# check missingness of covariates
sort(sapply(unlist(colnames_list$x), function(nm) {
	mean(is.na(df_unsplit[[nm]]))
}), decreasing=TRUE)

# select stimuli conditions
conditions_included = setdiff(df_unsplit[['stimuli']], 
'standard') # everything except 'standard'
df_unsplit = subset(df_unsplit, stimuli %in% conditions_included)
table(df_unsplit[['stimuli']])

# class label as stimuli
df_unsplit[['stimuli']] = ifelse(df_unsplit[['stimuli']]=='inattentive', 
1, 0) # inattentive is class 1
table(df_unsplit[['stimuli']])

# show response scale labels
levels(df_unsplit[['responseScale']])

# recode response scale
responsescale2group = setNames(c(3, 4, 5, 6, 7, 10, 10, 11, 11),
c('3 cat. (all labelled)', '4 cat. (all labelled)', '5 cat. (all labelled)', 
'6 cat. (all labelled)', '7 cat. (all labelled)', '10 cat. all labelled', 
'10 cat. end labelled', '11 cat. all labelled', '11 cat. end labelled'))
responsescale2group

# split by pointscale 
df_by_pointscale = with(new.env(), {
	pointscale_grouped = responsescale2group[df_unsplit[['responseScale']]]	
	pointscale_grouped = paste('ps', sprintf('%02d', pointscale_grouped), sep='')
	split(df_unsplit, f=pointscale_grouped)
})
lapply(df_by_pointscale, function(df) {
	table(df[['responseScale']], df[['stimuli']])
})

# check that splits make sense
lapply(df_by_pointscale, function(df) {
	t(sapply(df[, colnames_list$z], table))
})

# for 11-point scale, make it from 1 to 11
df_by_pointscale$ps11 = (function(df) {
	incremented = lapply(colnames(df), function(nm) {
		if(nm %in% colnames_list$z) {
			df[[nm]]+1
		} else {
			df[[nm]]
		} # if likert column, add 1
	})
	as.data.frame(setNames(incremented, colnames(df)))
})(df_by_pointscale$ps11)
t(sapply(df_by_pointscale$ps11[, colnames_list$z], table))

# per split, create list of xyz triple
processed_data = lapply(df_by_pointscale, function(df) {
	lapply(colnames_list, function(g) {
		df[, g]
	})
})
t(sapply(processed_data, names))

# save processed data
Sys.time()
saveRDS(processed_data, './data/analysis-preprocess.RDS')
devtools::session_info()
