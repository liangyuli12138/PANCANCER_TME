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
if(exists $ha{$a[2]}){for($i=3;$i<@a;$i++){$hb{$ha{$a[2]}}{$a[0]}{$i}+=$a[$i];}}};

for ($j=0;$j<@r;$j++){
for $i(sort {$a cmp $b} keys %{$hb{$r[$j]}}){
print "$r[$j]\t$r[$j]\t$i";
for $k(sort {$a <=> $b} keys %{$hb{$r[$j]}{$i}}){
$hb{$r[$j]}{$i}{$k}=$hb{$r[$j]}{$i}{$k}?$hb{$r[$j]}{$i}{$k}:0;
print "\t$hb{$r[$j]}{$i}{$k}"
}
;print "\n"
}
}

