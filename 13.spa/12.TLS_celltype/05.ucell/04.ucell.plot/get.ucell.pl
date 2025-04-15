open IN0,$ARGV[0];
<IN0>;
while(<IN0>){
chomp;
s/\"//g;
@a=split(/,/);
if($a[2]=~/Cluster/){$ha{$a[1]}=$a[2]}
}

open IN2,$ARGV[2];
$t=$ARGV[2];
open IN1,$ARGV[1];
<IN1>;
while(<IN1>){
s/\"//g;
@a=split(/,/);
if(exists $ha{$a[0]}){$hb{$ha{$a[0]}}+=$a[-2];$hc{$ha{$a[0]}}+=$a[-1];$hd{$ha{$a[0]}}++}
}

for $i(keys %hb){
$x=$hb{$i}/$hd{$i};$y=$hc{$i}/$hd{$i};
print "$t$i,$x,$y\n"
}
