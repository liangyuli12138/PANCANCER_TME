open IN0,$ARGV[0];
while(<IN0>){
chomp;@a=split;$ha{$a[0]}=$a[1]}

open IN1,$ARGV[1];
$t=<IN1>;print "$t";
while(<IN1>){
chomp;
$o=$_;
s/\"//g;
@a=split(/,/);
print "$a[0],$a[1],$a[2],$a[3],$a[4],$ha{$a[4]}\n" 
}
