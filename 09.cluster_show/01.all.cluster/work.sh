cut -d "," -f 1,16-21 pancancer.ref.0905.final.obs.csv > oridata/merge.all.all.level.sev.at

perl -e '$t=<>;chomp $t;@tt=split(/,/,$t);while(<>){chomp;@a=split/,/;for($i=15;$i<@a;$i++){$ha{$tt[$i]}{$a[$i]}=1}};for $i(keys %ha){print "$i";for $j(sort {$a cmp $b} key %{$ha{$i}}){print "\t$i"};print "\n"}' pancancer.ref.0905.final.obs.csv > pancancer.ref.0905.final.obs.csv.type

