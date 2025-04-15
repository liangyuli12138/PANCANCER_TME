open IN0,$ARGV[0];
while(<IN0>){chomp;@a=split;$ha{$a[0]}=$a[1]};

open IN1,$ARGV[1];<IN1>;
while(<IN1>){chomp;@a=split(/,/);$hb{$a[0]}="$a[-3],$ha{$a[-3]}"}

open IN2,$ARGV[2];$t=<IN2>;chomp $t;
print "$t,ori_celltype,target_celltype\n";
while(<IN2>){chomp;@a=split(/,/);print "$_,$hb{$a[0]}\n"}
