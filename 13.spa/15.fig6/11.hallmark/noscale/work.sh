perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;$ha{$_}=1};print "\"\"";open IN1,$ARGV[1];$t=<IN1>;chomp $t;$t=~s/\"//g;$t=~s/\_UCell//g;@x=split(/\,/,$t);for($i=1;$i<@x;$i++){if(exists $ha{$x[$i]}){$hb{$i}=1;print ",\"$x[$i]\""}};print "$x[$i]\n";while(<IN1>){chomp;s/\"//g;@a=split(/\,/);print "$a[1]";for($i=1;$i<@a;$i++){if(exists $hb{$i}){print ",$a[$i]"}};print "\n"}' hallmark.list hallmark.30.score.csv > hallmark.30.score.csv.filter

perl get.score.pl pancancer.icar.background.cell.obs hallmark.30.score.csv.filter hallmark.30.score.ori.csv > merge.gene.signature.score.csv &

