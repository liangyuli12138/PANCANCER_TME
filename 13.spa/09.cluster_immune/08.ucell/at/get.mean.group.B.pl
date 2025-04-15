open IN0,$ARGV[0];
while(<IN0>){
chomp;@a=split;
$t="at/$a[0].merge.at";
open IN,$t;
<IN>;
while(<IN>){
chomp;
@b=split(/,/);if(@b<9 || $b[3] !~/Cluster/){next};
$p=$a[0].$b[3];
#print "$p\n";
$ha{$p}++;
if($b[2] =~ /Lymphoid_B/ ){
$ha5{$p}+=$b[5];
$ha6{$p}+=$b[6];
$ha7{$p}+=$b[7];
$ha8{$p}+=$b[8];
}}}

print "cell,early,intermediate,Late,GO0002317\n";
open IN1,$ARGV[1];
<IN1>;
while(<IN1>){
chomp;@a=split(/,/);
if(exists $ha{$a[0]}){
$ha5{$a[0]}=$ha5{$a[0]}/$ha{$a[0]};
$ha6{$a[0]}=$ha6{$a[0]}/$ha{$a[0]};
$ha7{$a[0]}=$ha7{$a[0]}/$ha{$a[0]};
$ha8{$a[0]}=$ha8{$a[0]}/$ha{$a[0]};
print "$a[0],$ha5{$a[0]},$ha6{$a[0]},$ha7{$a[0]},$ha8{$a[0]}\n"
}
else{print "$a[0],0,0,0,0\n"}
}
