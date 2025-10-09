# real analysis data
wget https://osf.io/download/kv9xj/ -O auxiliary-data-wide.csv
wget https://osf.io/download/hdmvs/ -O "experiments1&4-results-to-analysis.RDS"

# HSQ data
wget http://openpsychometrics.org/_rawdata/HSQ.zip
unzip -j HSQ.zip HSQ/data.csv
mv data.csv openpsychometrics-hsq.csv # rename
rm HSQ.zip
