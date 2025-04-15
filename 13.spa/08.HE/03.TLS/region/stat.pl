open IN0,$ARGV[0];
while(<IN0>){chomp;@a=split;$ha{$a[0]}=$a[1]};

open IN1,$ARGV[1];
while(<IN1>){
chomp;
@a=split(/,/);
$hb{$a[0]}=$a[1];
}

#$s=$ARGV[3];

open IN2,$ARGV[2];
<IN2>;
while(<IN2>){
#s/\"//g;
chomp;
@a=split(/,/);
if(exists $hb{$a[0]} && exists $ha{$a[-3]}){$hc{$hb{$a[0]}}{$ha{$a[-3]}}++;$hd{$hb{$a[0]}}{$a[-3]}++};
#if(!exists $ha{$a[-3]}){print "$a[0]\t$a[-3]\n";}
}

#for $i(sort {$a cmp $b} keys %hc){for $j(sort {$a cmp $b} keys %{$hc{$i}}){print "$i\t$j\t$hc{$i}{$j}\t$j\n"}}

for $i(sort {$a cmp $b} keys %hd){
for $j(sort {$a cmp $b} keys %{$hd{$i}})
{
print "$i\t$j\t$hd{$i}{$j}\n"
}
}
