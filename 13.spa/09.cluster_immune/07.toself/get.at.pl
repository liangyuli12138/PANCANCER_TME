open IN0,$ARGV[0];
$t=$ARGV[0];
$id=$1 if($t=~/at\/(\S+)\.ex/);
<IN0>;
while(<IN0>){
chomp;
if(/other/){next};
s/\"//g;
@a=split(/,/);
$a[2]=~s/Lymphoid_B_naive/Lymphoid_B/;$a[2]=~s/Lymphoid_B_memory/Lymphoid_B/;
$ha{$a[3]}.="$a[1]\n";
$hb{$a[3]}.="$a[1],$a[2]\n";
}

for $i(keys %ha){
$ip=$id."_".$i;
open OUT0,">at/$ip.input";
open OUT1,">at/$ip.at";
print OUT0 "cell\n$ha{$i}";
print OUT1 "cell,tls_type\n$hb{$i}";
}
