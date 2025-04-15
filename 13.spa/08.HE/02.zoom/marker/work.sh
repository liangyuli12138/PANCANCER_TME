perl -e 'print "\"\",\"id\",\"region\"\n";<>;while(<>){chomp;$n++;@a=split(/,/);if(/immune/){print "\"$n\",\"$a[0]\",\"immune\"\n"}else{print "\"$n\",\"$a[0]\",\"others\"\n"}}' D01972D1_region.csv > D01972D1_region.csv.at

perl -e 'print "cell,region\n";<>;while(<>){chomp;$n++;@a=split(/,/);if(/immune/){print "$a[0],immune\n"}else{print "$a[0],others\n"}}' D01972D1_region.csv > D01972D1_region.csv.at

