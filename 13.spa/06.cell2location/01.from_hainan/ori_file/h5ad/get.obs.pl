open IN2,$ARGV[2];
while(<IN2>){
chomp;
@a=split;
$hm{$a[0]}=$a[1];
}
open IN0,$ARGV[0];
<IN0>;
while(<IN0>){
chomp;
s/\"//g;
@a=split(/,/);
$ha{$a[1]}="$a[2],$hm{$a[2]},$a[3]"
}
open IN1,$ARGV[1];
$t=<IN1>;
chomp $t;
print "$t,celltype,celltype_merge,max_score\n";
while(<IN1>){
chomp;
@a=split(/,/);
print "$_,$ha{$a[0]}\n"
}
