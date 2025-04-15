perl -e 'open IN0,$ARGV[0];<IN0>;while(<IN0>){chomp;s/\t+$//;@a=split;for($i=1;$i<@a;$i++){$ha{$a[$i]}=$a[0]}};open IN1,$ARGV[1];while(<IN1>){chomp;open IN,$_;$id=$1 if(/out\/(\S+)\//);$t=<IN>;chomp $t;print "tissues,tls_marker,$t\n";while(<IN>){chomp;@a=split(/,/);if(exists $ha{$a[1]}){print "$a[0],$id,$ha{$a[1]},$_\n"}}}' tls.marker.list immune.diff.list > tls.marker.diff.csv

