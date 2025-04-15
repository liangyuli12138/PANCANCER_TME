perl stat.pl cell.list region.list sn.list > all.he.regeion.cell.stat.xls

perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;@a=split;$ha{$a[0]}{$a[1]}=$a[2]};open IN1,$ARGV[1];while(<IN1>){chomp;@a=split;$input="region.old/$a[0]"."_region.csv";open IN,$input;$out="region/$a[0]"."_region.csv";open OUT,">$out";while(<IN>){chomp;@b=split(/,/);if(exists $ha{$a[0]}{$b[1]}){print OUT "$b[0],$ha{$a[0]}{$b[1]},$b[2],$b[3]\n"}else{if($a[0] eq "B02324E3" && $b[1] eq "fat"){next};print OUT "$_\n"}}}' change.list sn.list


perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;@a=split;$ha{$a[0]}=$a[1]};open IN1,$ARGV[1];while(<IN1>){chomp;@a=split(/\//);open IN,$_;open OUT,">region_merge/$a[-1]";$t=<IN>;print OUT "$t";while(<IN>){chomp;@b=split(/,/);if(exists $ha{$b[1]}){print OUT "$b[0],$ha{$b[1]},$b[2],$b[3]\n"}}}' merge.list region.list



