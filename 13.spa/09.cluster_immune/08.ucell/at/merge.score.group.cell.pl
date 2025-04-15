open IN0,$ARGV[0];
<IN0>;
while(<IN0>){chomp;@a=split(/,/);$ha{$a[0]}="$a[4],$a[5],$a[6],$a[7]"}

open IN1,$ARGV[1];
<IN1>;
while(<IN1>){chomp;@a=split(/,/);$hb{$a[1]}=$a[2]}

open IN2,$ARGV[2];
$t=<IN2>;chomp $t;
print "$t,\"cellgroup\",\"early_UCell\",\"intermediate_UCell\",\"Late_UCell\",\"GO0002317_UCell\"\n";
while(<IN2>){chomp;@a=split(/,/);print "$_,$hb{$a[1]},$ha{$a[1]}\n"}
