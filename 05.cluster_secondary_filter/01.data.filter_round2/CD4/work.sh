perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;$ha{$_}=1};open IN1,$ARGV[1];while(<IN1>){chomp;if(exists $ha{$_}){}else{print "$_\n"}}' monify1.cd4.cluster8.list filter.CD8.list.input > filter.CD8.list.monify2.input

