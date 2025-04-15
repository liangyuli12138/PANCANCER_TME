open IN0,$ARGV[0];
<IN0>;
while(<IN0>){
chomp;
s/\"//g;
@a=split(/,/);
if($a[4]=~/ICM/){$ha{$a[0]}{$a[4]}=1}
if($a[6] ne "other"){$ha{$a[0]}{$a[6]}=1}
};

open IN1,$ARGV[1];
<IN1>;
while(<IN1>){
if(/^id/){next};
chomp;$c=$_;
s/\"//g;
@a=split(/,/);
for $i(keys %ha){
for $j(keys %{$ha{$i}}){
@x=split(/\_/,$a[0]);
@y=split(/\_/,$i);
$l=int(sqrt(($x[0]-$y[0])*($x[0]-$y[0])+($x[1]-$y[1])*($x[1]-$y[1])));
$hb{$a[0]}{$j}{$l}=$i;
$hc{$a[0]}=$c;
}}}

for $i(keys %hb){
for $j(keys %{$hb{$i}}){
$n=0;$m=0;$t="";
for $k(sort {$a<=>$b} keys %{$hb{$i}{$j}}){
$n++;if($n<=5){$m+=$k;$t.="$hb{$i}{$j}{$k}|"}else{$n=$n-1;last}
}
$m=$m/$n;
print "$i,$j,$m,$n,$t,$hc{$i}\n"
}}
