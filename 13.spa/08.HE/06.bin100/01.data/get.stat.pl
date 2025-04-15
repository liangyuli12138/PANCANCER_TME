open IN0,$ARGV[0];
while(<IN0>){
chomp;
$id=$1 if(/data\/(\S+)\_cellbin/);
open IN,$_;
<IN>;
while(<IN>){
if(/Unknown/ || /EC_Glomerular/ || /EC_Sinusoidal/ || /Hepatocyte/){next};
@a=split(/,/);
$x=int($a[-5]/200)*200;$y=int($a[-4]/200)*200;
$z="$id"."_"."$x"."_"."$y";
$ha{$z}{$a[-3]}++;$haa{$z}++;
$hb{$z}{$a[-2]}++;$hbb{$z}++;
$hc{$a[-3]}=1;$hd{$a[-2]}=1;
}}

open OUT0,">$ARGV[1]";
for $i(sort {$a cmp $b} keys %hc){
print OUT0 ",$i"
}
print OUT0 "\n";
for $j(sort {$a cmp $b} keys %ha){
print OUT0 "$j";
for $i(sort {$a cmp $b} keys %hc){
$o=$ha{$j}{$i}?$ha{$j}{$i}/$haa{$j}:0;
print OUT0 ",$o";
}
print OUT0 "\n";
}

open OUT0,">$ARGV[2]";
print OUT0 ",total";
for $i(sort {$a cmp $b} keys %hc){
print OUT0 ",$i"
}
print OUT0 "\n";
for $j(sort {$a cmp $b} keys %ha){
print OUT0 "$j,$haa{$j}";
for $i(sort {$a cmp $b} keys %hc){
$o=$ha{$j}{$i}?$ha{$j}{$i}:0;
print OUT0 ",$o";
}
print OUT0 "\n";
}

open OUT0,">$ARGV[3]";
for $i(sort {$a cmp $b} keys %hd){
print OUT0 ",$i"
}
print OUT0 "\n";
for $j(sort {$a cmp $b} keys %hb){
print OUT0 "$j";
for $i(sort {$a cmp $b} keys %hd){
$o=$hb{$j}{$i}?$hb{$j}{$i}/$hbb{$j}:0;
print OUT0 ",$o";
}
print OUT0 "\n";
}

open OUT0,">$ARGV[4]";
print OUT0 ",total";
for $i(sort {$a cmp $b} keys %hd){
print OUT0 ",$i"
}
print OUT0 "\n";
for $j(sort {$a cmp $b} keys %hb){
print OUT0 "$j,$hbb{$j}";
for $i(sort {$a cmp $b} keys %hd){
$o=$hb{$j}{$i}?$hb{$j}{$i}:0;
print OUT0 ",$o";
}
print OUT0 "\n";
}

