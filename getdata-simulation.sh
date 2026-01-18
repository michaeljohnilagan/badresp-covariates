# HSQ data
wget http://openpsychometrics.org/_rawdata/HSQ.zip
unzip -j HSQ.zip HSQ/data.csv
mv data.csv ./data/openpsychometrics-hsq.csv # rename
rm HSQ.zip
