open IN0,$ARGV[0];
<IN0>;
while(<IN0>){
chomp;
s/\"//g;
@a=split(/,/);
if($a[2]=~/Lymphoid_B/){$ha{$a[1]}="$a[2],$a[3]"}
}

open IN2,$ARGV[2];
$t=$ARGV[2];
open IN1,$ARGV[1];
<IN1>;
while(<IN1>){
s/\"//g;chomp;
@a=split(/,/);
if(exists $ha{$a[0]}){print "$t,$_,$ha{$a[0]}\n"}
}

