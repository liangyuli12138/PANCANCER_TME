perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;@a=split;$ha{$a[0]}{$a[1]}=$a[3]};open IN1,$ARGV[1];while(<IN1>){chomp;open IN,$_;$id=$1 if(/at\/(\S+)\.immune/);<IN>;while(<IN>){chomp;s/\"//g;@a=split(/,/);if($a[2]=~/Clu/){print "$id\t$a[1]\t$a[2]\t$id$a[2]\t$ha{$id}{$a[1]}\n"}}}' merge.stat.csv immune.list > merge.stat.csv.link

perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;@a=split;$ha{$a[0]}=1};open IN1,$ARGV[1];while(<IN1>){chomp;@a=split(/\t/);if(exists $ha{$a[3]}){}else{print "$_\n"}}' filter.man.cluster.list merge.stat.csv.link > merge.stat.csv.link.filter

