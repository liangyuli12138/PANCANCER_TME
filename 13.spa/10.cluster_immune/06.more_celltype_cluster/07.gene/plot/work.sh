perl -e 'while(<>){chomp;@a=split(/,/);if(exists $ha{$a[1]}){next};$ha{$a[1]}=1;s/\d//;print "$_\n"}' merge.csv > merge.csv.filter

