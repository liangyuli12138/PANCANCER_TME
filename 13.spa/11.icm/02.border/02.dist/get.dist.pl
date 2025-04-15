<>;
while(<>){
chomp;
s/\"//g;
@a=split(/,/);
$ha{$a[1]}=$a[2]};

for $i(keys %ha){for $j(keys %ha){
if($i eq $j || $ha{$j} ne "other"){next};
@x=split(/\_/,$i);
@y=split(/\_/,$j);
$l=int(sqrt(($x[0]-$y[0])*($x[0]-$y[0])+($x[1]-$y[1])*($x[1]-$y[1])));
if(!exists $hx{$i}){$hx{$i}=$l}else{if($l<$hx{$i}){$hx{$i}=$l}}
}}

for $i(keys %ha){for $j(keys %ha){
if($i eq $j || $ha{$j} ne "malignant"){next};
@x=split(/\_/,$i);
@y=split(/\_/,$j);
$l=int(sqrt(($x[0]-$y[0])*($x[0]-$y[0])+($x[1]-$y[1])*($x[1]-$y[1])));
if(!exists $hy{$i}){$hy{$i}=$l}else{if($l<$hy{$i}){$hy{$i}=$l}}
}}

for $i(keys %hx){print "$i\t$ha{$i}\t$hx{$i}\t$hy{$i}\n"}
