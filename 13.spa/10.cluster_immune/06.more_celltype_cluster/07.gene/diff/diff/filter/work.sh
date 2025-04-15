ls diff.*filter |perl -e 'while(<>){chomp;$id=$1 if(/diff\.(\S+)\.csv/);open IN,$_;while(<IN>){@a=split(/,/);print "\"$id\",\"$a[1]\"\n"}}' > merge.csv

