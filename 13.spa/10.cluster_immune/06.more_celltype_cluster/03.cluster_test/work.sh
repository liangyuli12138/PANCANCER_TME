perl -e 'for($i=3;$i<=20;$i++){$a=`cat temp.py`;$a=~s/aaaa/$i/g;print "$a"}' >> cluster.py 
perl -e 'for($i=3;$i<=20;$i++){$p="out/immune.cluster.$i.r1.5.obs";open IN,$p;<IN>;while(<IN>){@a=split(/,/);$ha{$a[0]}{$i}=$a[1]}};for $j(keys %ha){print "$j";for($i=3;$i<=20;$i++){print "\t$ha{$j}{$i}"};print "\n"}' > merge.all.test.obs.xls

