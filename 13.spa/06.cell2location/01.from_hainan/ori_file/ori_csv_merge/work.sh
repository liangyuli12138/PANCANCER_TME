perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;@a=split;$ha{$a[0]}=$a[1]};open IN1,$ARGV[1];while(<IN1>){chomp;@b=split(/\//);$o=$b[2];open OUT,">$o";open IN,$_;$t=<IN>;print OUT "$t";while(<IN>){chomp;@a=split(/,/);$a[2]=~s/\"//g;s/$a[2]/$ha{$a[2]}/;print OUT "$_\n"}}' cluster.list ori.list

