# set up reproducible environment
if [ ! -d './renv/library' ]; then
R -e 'renv::restore()'
fi

# not simulation, not analysis
if [ ! -f './tests.Rout' ]; then 
R --no-save --no-restore --file=tests/tests.R > tests.Rout 2>&1
fi
if [ ! -f './figs-other.Rout' ]; then 
R --no-save --no-restore --file=scripts/figs-other.R > figs-other.Rout 2>&1
fi

# get data for simulation
if [ ! -f './data/openpsychometrics-hsq.csv' ]; then
bash scripts/simulation/01-getdata-simulation.sh
fi

# run simulation
if [ ! -f 'simulation.RData' ]; then
R --no-save --no-restore --file=scripts/simulation/02-simulation.R > 02-simulation.Rout 2>&1
R --no-save --no-restore --file=scripts/simulation/03-figs-simulation.R > 03-figs-simulation.Rout 2>&1
fi

# get data for analysis
if [ ! -f './data/auxiliary-data-wide.csv' ] || [ ! -f './data/experiments1&4-results-to-analysis.RDS' ] ; then
bash scripts/analysis/01-getdata-analysis.sh
fi

# run analysis
if [ ! -f 'analysis.RData' ]; then
R --no-save --no-restore --file=scripts/analysis/02-analysis-preprocess.R > 02-analysis-preprocess.Rout 2>&1
R --no-save --no-restore --file=scripts/analysis/03-analysis-main.R > 03-analysis-main.Rout 2>&1
R --no-save --no-restore --file=scripts/analysis/04-figs-analysis.R > 04-figs-analysis.Rout 2>&1
fi
