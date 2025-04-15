@r=("malignant","epithelium","immune","fibroblast","vessel","muscle","fat");
open IN0,$ARGV[0];
while(<IN0>){
chomp;@a=split;$ha{$a[0]}=$a[1]};

open IN1,$ARGV[1];
$t=<IN1>;
chomp $t;@t=split(/\t/,$t);print "$t\n";

while(<IN1>){
chomp;
s/\t+$//;
@a=split;
if(exists $ha{$a[2]}){for($i=3;$i<@a;$i++){$hb{$a[0]}{$ha{$a[2]}}{$i}+=$a[$i];}}};

for $i(sort {$a cmp $b} keys %hb)
{
for ($j=0;$j<@r;$j++){
print "$i\t$i\t$r[$j]";
for $k(sort {$a <=> $b} keys %{$hb{$i}{$r[$j]}}){
$hb{$i}{$r[$j]}{$k}=$hb{$i}{$r[$j]}{$k}?$hb{$i}{$r[$j]}{$k}:0;
print "\t$hb{$i}{$r[$j]}{$k}"};print "\n"}}
