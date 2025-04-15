open IN0,$ARGV[0];
<IN0>;
while(<IN0>){
chomp;
s/\"//g;
@a=split(/,/);
$ha{$a[1]}=$a[2]};

open IN1,$ARGV[1];
<IN1>;
while(<IN1>){
chomp;
s/\"//g;
@a=split(/,/);
$hb{$a[1]}=$a[2]};

for $i(keys %ha){for $j(keys %hb){
if($i eq $j || $hb{$j} ne "other"){next};
@x=split(/\_/,$i);
@y=split(/\_/,$j);
$l=int(sqrt(($x[0]-$y[0])*($x[0]-$y[0])+($x[1]-$y[1])*($x[1]-$y[1])));
if(!exists $hx{$i}){$hx{$i}=$l}else{if($l<$hx{$i}){$hx{$i}=$l}}
}}

for $i(keys %ha){for $j(keys %hb){
if($i eq $j || $hb{$j} ne "malignant"){next};
@x=split(/\_/,$i);
@y=split(/\_/,$j);
$l=int(sqrt(($x[0]-$y[0])*($x[0]-$y[0])+($x[1]-$y[1])*($x[1]-$y[1])));
if(!exists $hy{$i}){$hy{$i}=$l}else{if($l<$hy{$i}){$hy{$i}=$l}}
}}

for $i(keys %hx){print "$i\t$ha{$i}\t$hx{$i}\t$hy{$i}\n"}
