open IN0,$ARGV[0];
<IN0>;
while(<IN0>){
chomp;
@a=split(/,/);
if($a[3]>2 && $a[5]<0.05 && $a[6]>0.05){$ha{$a[3]}=$_}
}

for $i(sort {$b <=> $a} keys %ha){if($i>1){print "$ha{$i}\n"}}
