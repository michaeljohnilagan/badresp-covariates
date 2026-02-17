# set up reproducible environment
date; R -e 'renv::restore()'; date

# not simulation, not analysis
date; R --no-save --no-restore --file=scripts/tests.R > tests.Rout 2>&1; date
date; R --no-save --no-restore --file=scripts/figs-other.R > figs-other.Rout 2>&1; date

# simulation
bash scripts/simulation/01-getdata-simulation.sh
date; R --no-save --no-restore --file=scripts/simulation/02-simulation.R > 02-simulation.Rout 2>&1; date
R --no-save --no-restore --file=scripts/simulation/03-figs-simulation.R > 03-figs-simulation.Rout 2>&1

# analysis
bash scripts/analysis/01-getdata-analysis.sh
date; R --no-save --no-restore --file=scripts/analysis/02-analysis-preprocess.R > 02-analysis-preprocess.Rout 2>&1; date
date; R --no-save --no-restore --file=scripts/analysis/03-analysis-main.R > 03-analysis-main.Rout 2>&1; date
R --no-save --no-restore --file=scripts/analysis/04-figs-analysis.R > 04-figs-analysis.Rout 2>&1
