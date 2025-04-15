open IN0,$ARGV[0];
<IN0>;
while(<IN0>){
chomp;s/\"//g;
@a=split(/,/);
$hz{$a[1]}=1;
$hn{$a[1]}++;
$hm{$a[1]}+=$a[5];
$ha{$a[1]}{$a[2]}++;
$hb{$a[1]}{$a[2]}+=$a[5];
}

open IN1,$ARGV[1];
<IN1>;
while(<IN1>){
chomp;s/\"//g;
@a=split(/,/);
$hz{$a[1]}=1;
$hnn{$a[1]}++;
$hmm{$a[1]}+=$a[5];
$haa{$a[1]}{$a[2]}++;
$hbb{$a[1]}{$a[2]}+=$a[5];
}

open OUT0,">$ARGV[2]";
open OUT1,">$ARGV[3]";

print OUT0 "cell";
print OUT1 "cell";

for $i(sort {$a cmp $b} keys %hz){print OUT0 "\t$i";print OUT1 "\t$i"};print OUT0 "\n";print OUT1 "\n";

for $j(sort {$a cmp $b} keys %hz){

print OUT0 "$j";print OUT1 "$j";

for $i(sort {$a cmp $b} keys %hz){


$x=sprintf "%.2f",($ha{$j}{$i}/$hn{$j})/($haa{$j}{$i}/$hnn{$j});
$y=sprintf "%.2f",($hb{$j}{$i}/$hm{$j})/($hbb{$j}{$i}/$hmm{$j});

#$x=$ha{$j}{$i}/$hn{$j};$y=$hb{$j}{$i}/$hm{$j};

print OUT0 "\t$x";
print OUT1 "\t$y";

}
print OUT0 "\n";
print OUT1 "\n";
}
