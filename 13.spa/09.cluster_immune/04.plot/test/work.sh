perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;s/\"//g;@a=split(/,/);if($a[2]=~/group/){$ha{$a[1]}=$a[2]}};open IN1,$ARGV[1];while(<IN1>){chomp;s/\"//g;@a=split(/,/);if(exists $ha{$a[1]}){print "$_,$ha{$a[1]}\n"}}' D01972D1.group.at D01972D1.ori.at > D01972D1.all.immune.csv

