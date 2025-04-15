#perl -e 'while(<>){s/\t+/\t/g;s/\t$//;print "$_"}' marker.list.ori > marker.list
perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;$ha{$_}=1};open IN1,$ARGV[1];while(<IN1>){chomp;@a=split;print "$a[0]";for($i=1;$i<@a;$i++){if(exists $ha{$a[$i]}){print "\t$a[$i]"}};print "\n"}' pancancer.final.0314.var.csv.input.add marker.list.ori > marker.list


perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;@a=split;$n++;$ha{$a[0]}=$n};open IN1,$ARGV[1];while(<IN1>){chomp;@a=split;s/[^\t]+\t//;$hb{$a[0]}=$_};open IN2,$ARGV[2];while(<IN2>){chomp;@a=split(/,/);$hc{$a[-1]}{$a[0]}=1};for $i(keys %hc){for $j(sort {$ha{$a} <=> $ha{$b}} %ha){if(exists $hc{$i}{$j}){print "$i\t$j\t$hb{$j}\n"}}}'  cell.groups_pri.list marker.list cell.groups_sev.list > cluster.marker.lsit

perl get.py.pl cell.groups_pri.list marker.list cell.groups_sev.list > plot.py

perl -e 'while(<>){chomp;if(/\[(.+)\],$/){$a=$1;print "$a,\n"}}' plot.py |les 

