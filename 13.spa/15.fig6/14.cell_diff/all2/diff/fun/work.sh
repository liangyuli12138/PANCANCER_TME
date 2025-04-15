perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;$ha{$_}=1};open IN1,$ARGV[1];while(<IN1>){chomp;$id=$1 if(/diff\.(\S+)\_Lym/);if(exists $ha{$id}){open IN,$_;$n=0;<IN>;while(<IN>){chomp;@a=split(/,/);$n++;if(/,RPS/ || /,RPL/){next};if($n<=100){$hb{$a[1]}.="$id"."|";$hc{$a[1]}++;}}}};for $i(keys %hb){print "$i\tAll\t$hc{$i}\t$hb{$i}\n"}' pick.celltype.list diff.list |sort -k 3,3nr > all.diff.gene.list

perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;$ha{$_}=1};open IN1,$ARGV[1];while(<IN1>){chomp;$id=$1 if(/diff\.(\S+)\_Lym/);if(exists $ha{$id}){open IN,$_;$n=0;<IN>;while(<IN>){chomp;@a=split(/,/);$n++;if(/,RPS/ || /,RPL/){next};if($n<=100){print "$a[1]\t$id\n"}}}}' pick.celltype.list diff.list > all.for.enrich.gene.list2

perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;$ha{$_}=1};open IN1,$ARGV[1];while(<IN1>){chomp;$id=$1 if(/diff\.(\S+)\_Lym/);if(exists $ha{$id}){open IN,$_;$n=0;<IN>;while(<IN>){chomp;@a=split(/,/);$n++;if(/,RPS/ || /,RPL/){next};if($n<=100){print "$a[1]\t$id\n"}}}}' pick.celltype.list diff.list2 > all.for.enrich.gene.list3

