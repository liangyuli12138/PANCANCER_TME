perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;@a=split;$ha{$a[0]}=$a[1]};open IN1,$ARGV[1];while(<IN1>){chomp;@a=split;for($i=2;$i<@a;$i++){$hb{$a[$i]}.="$a[1]|";if(exists $ha{$a[1]}){$hc{$a[$i]}.="$ha{$a[1]}|"};}};for($i=1;$i<=3;$i++){$f="icm$i.diff.csv";open IN,$f;<IN>;$n=0;while(<IN>){chomp;@a=split(/,/);if($n<2000 && $a[3]>1 && $a[5]<0.05){$x="ICM$i";if(exists $hc{$a[1]} && $hc{$a[1]}=~/$x/){print "ICM$i,$hb{$a[1]},$hc{$a[1]},$_\n";$n++}}}}' icm.list cluster.marker.lsit |les

perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;@a=split;$ha{$a[0]}=$a[1]};open IN1,$ARGV[1];while(<IN1>){chomp;@a=split;for($i=2;$i<@a;$i++){$hb{$a[$i]}.="$a[1]|";if(exists $ha{$a[1]}){$hc{$a[$i]}.="$ha{$a[1]}|"};}};for($i=1;$i<=3;$i++){$f="icm$i.diff.csv";open IN,$f;<IN>;$n=0;while(<IN>){chomp;@a=split(/,/);if($n<100 && $a[3]>1 && $a[5]<0.05){$x="ICM$i";print "ICM$i,$hb{$a[1]},$hc{$a[1]},$_\n";$n++}}}' icm.list cluster.marker.lsit > cluster.marker.lsit.forplot.in

perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;@a=split;$ha{$a[0]}=$a[1]};open IN1,$ARGV[1];while(<IN1>){chomp;@a=split;for($i=2;$i<@a;$i++){$hb{$a[$i]}.="$a[1]|";if(exists $ha{$a[1]}){$hc{$a[$i]}.="$ha{$a[1]}|"};}};for($i=1;$i<=3;$i++){$f="icm$i.diff.csv";open IN,$f;<IN>;$n=0;while(<IN>){chomp;@a=split(/,/);if($n<100 && $a[3]>1 && $a[5]<0.05){$x="ICM$i";print "ICM$i,$hb{$a[1]},$hc{$a[1]},$_\n";$n++}}}' icm.list cluster.marker.lsit > cluster.marker.lsit.forplot

awk -F "," '$3~/ICM/ {print $1","$5}' cluster.marker.lsit.forplot > cluster.marker.lsit.forplot.input

