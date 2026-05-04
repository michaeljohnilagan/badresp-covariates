#!/bin/bash

#SBATCH --time=40:00:00
#SBATCH --mem_per-cpu=500M
#SBATCH --cpus-per-task=1
#SBATCH --job-name='milagan_badresp-covariates'

# set up reproducible environment
echo 'setting up reproducible environment...'
if [ ! -d './renv/library' ]; then
R -e 'renv::restore()'
fi

# not simulation, not analysis
echo 'running tests...'
if [ ! -f './tests.Rout' ]; then 
R --no-save --no-restore --file=tests/tests.R > tests.Rout 2>&1
fi
if [ ! -f './figs-other.Rout' ]; then 
R --no-save --no-restore --file=scripts/figs-other.R > figs-other.Rout 2>&1
fi

# get data for simulation
echo 'retrieving data for simulation study...'
if [ ! -f './data/openpsychometrics-hsq.csv' ]; then
bash scripts/simulation/01-getdata-simulation.sh
fi

# run simulation
echo 'running simulation study...'
if [ ! -f 'simulation.RData' ]; then
R --no-save --no-restore --file=scripts/simulation/02-simulation.R > 02-simulation.Rout 2>&1
R --no-save --no-restore --file=scripts/simulation/03-figs-simulation.R > 03-figs-simulation.Rout 2>&1
fi

# get data for analysis
echo 'retrieving data for real data analysis...'
if [ ! -f './data/auxiliary-data-wide.csv' ] || [ ! -f './data/experiments1&4-results-to-analysis.RDS' ] ; then
bash scripts/analysis/01-getdata-analysis.sh
fi

# run analysis
echo 'running real data analysis...'
if [ ! -f 'analysis.RData' ]; then
R --no-save --no-restore --file=scripts/analysis/02-analysis-preprocess.R > 02-analysis-preprocess.Rout 2>&1
R --no-save --no-restore --file=scripts/analysis/03-analysis-main.R > 03-analysis-main.Rout 2>&1
R --no-save --no-restore --file=scripts/analysis/04-figs-analysis.R > 04-figs-analysis.Rout 2>&1
fi
