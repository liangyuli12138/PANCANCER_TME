cat out/* > all.out.merge &
perl -e 'while(<>){chomp;$id=$1 if(/out\/(\S+)\.all/);open IN,$_;while(<IN>){chomp;@a=split(/,/);if($a[1]=~/ICM/ && $a[-3]=~/ICM/){print "$id,$_\n"}}}' out.list > merge.icm.out.csv &

perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;@a=split;$ha{$a[0]}=$a[1]};print "tissue,sample,id,icmv1,icmv2,distance\n";open IN1,$ARGV[1];while(<IN1>){chomp;@a=split(/,/);print "$ha{$a[0]},$a[0],$a[1],$a[2],$a[-3],$a[3]\n"}' all.sn.list merge.icm.out.csv > merge.icm.out.icm.csv &

perl -e 'while(<>){chomp;@a=split(/,/);$ha{$a[3]}{$a[4]}++;$hb{$a[3]}{$a[4]}+=$a[5]};for $i(keys %ha){for $j(keys %{$ha{$i}}){$o=$hb{$i}{$j}/$ha{$i}{$j};print "$i\t$j\t$o\n"}}' merge.icm.out.icm.csv > merge.icm.out.icm.csv.stat &

