open IN0,$ARGV[0];
<IN0>;
while(<IN0>){
chomp;
@a=split(/,/);
if($a[-2]=~/Lymphoid/){$ha{$a[0]}=$a[-2]}}

open IN1,$ARGV[1];
<IN1>;
while(<IN1>){
chomp;s/\"//g;
@a=split(/,/);
@b=split(/\./,$a[0]);
if(exists $ha{$b[-1]}){$hb{$b[-1]}.= "$ha{$b[-1]},$_"}
}

open IN2,$ARGV[2];
<IN2>;
while(<IN2>){
chomp;s/\"//g;
@a=split(/,/);
@b=split(/\./,$a[0]);
if(exists $ha{$b[-1]}){$hb{$b[-1]}.= ",$a[1],$a[2],$a[3],$a[4]"}
}

for $i(sort {$a cmp $b} keys %hb){
print "$hb{$i}\n"
}
