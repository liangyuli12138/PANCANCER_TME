perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;@a=split;$ha{$a[0]}=$a[1]};open IN1,$ARGV[1];$t=<IN1>;print "$t";while(<IN1>){chomp;@a=split(/,/);print "$_","_","$ha{$a[1]}\n"}' type.list ../../celltype_merge_zz/all.cell.dist.csv.type > all.cell.dist.csv.type
perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;@a=split;$ha{$a[0]}=$a[1]};open IN1,$ARGV[1];$t=<IN1>;print "$t";while(<IN1>){chomp;@a=split(/,/);print "$_","_","$ha{$a[1]}\n"}' type.list ../../celltype_merge_zz/all.cell.dist.csv.type.ICM > all.cell.dist.csv.type.ICM

perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;@a=split;$ha{$a[0]}=$a[1]};open IN1,$ARGV[1];$t=<IN1>;while(<IN1>){chomp;@a=split(/,/);$hb{$ha{$a[1]}}.="$_"."_"."$ha{$a[1]}\n"};for $i(keys %hb){open OUT,">all.cell.dist.csv.type.$i";print OUT "$t","$hb{$i}"}' type.list ../../celltype_merge_zz/all.cell.dist.csv.type &
perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;@a=split;$ha{$a[0]}=$a[1]};open IN1,$ARGV[1];$t=<IN1>;while(<IN1>){chomp;@a=split(/,/);$hb{$ha{$a[1]}}.="$_"."_"."$ha{$a[1]}\n"};for $i(keys %hb){open OUT,">all.cell.dist.csv.type.ICM.$i";print OUT "$t","$hb{$i}"}' type.list ../../celltype_merge_zz/all.cell.dist.csv.type.ICM &


