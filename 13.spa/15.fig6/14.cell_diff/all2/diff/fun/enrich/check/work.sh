perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;@a=split(/,/);$ha{$a[0]}=1};open IN1,$ARGV[1];while(<IN1>){chomp;@a=split;if(exists $ha{$a[0]}){print "$_\n"}}' Proteoglycans.gene all.diff.gene.list > Proteoglycans.gene.cell

