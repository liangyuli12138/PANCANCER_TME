open IN0,$ARGV[0];
<IN0>;
while(<IN0>){
chomp;
s/\"//g;
@a=split(/,/);
$ha{$a[0]}=$a[4]};

open IN1,$ARGV[1];
<IN1>;
while(<IN1>){
chomp;
s/\"//g;
@a=split(/,/);
$hb{$a[0]}=$a[4]};

for $i(keys %ha){for $j(keys %hb){
if($i eq $j){next};
if($hb{$j} eq "other"){
@x=split(/\_/,$i);
@y=split(/\_/,$j);
$l=int(sqrt(($x[0]-$y[0])*($x[0]-$y[0])+($x[1]-$y[1])*($x[1]-$y[1])));
if(!exists $hx{$i}){$hx{$i}=$l}else{if($l<$hx{$i}){$hx{$i}=$l}}
}
elsif($hb{$j} eq "malignant"){
@x=split(/\_/,$i);
@y=split(/\_/,$j);
$l=int(sqrt(($x[0]-$y[0])*($x[0]-$y[0])+($x[1]-$y[1])*($x[1]-$y[1])));
if(!exists $hy{$i}){$hy{$i}=$l}else{if($l<$hy{$i}){$hy{$i}=$l}}
}
}}

for $i(keys %hx){print "$i\t$ha{$i}\t$hx{$i}\t$hy{$i}\n"}
