perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;@a=split(/,/);$ha{$a[1]}=$a[0]};open IN1,$ARGV[1];while(<IN1>){chomp;@a=split(/,/);if(exists $ha{$a[1]}){print "$ha{$a[1]},$_\n"}}' ref.gene.list.input ../diff.Lymphoid1.csv > diff.Lymphoid1.csv.stat

