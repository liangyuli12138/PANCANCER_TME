while(<>){
if(/immue_region/){next};
chomp;
@a=split(/,/);
if(!exists $hxa{$a[1]}){$hxa{$a[1]}=$a[2];$hxb{$a[1]}=$a[2]}else{$hxa{$a[1]}=$hxa{$a[1]}>$a[2]?$a[2]:$hxa{$a[1]};$hxb{$a[1]}=$hxb{$a[1]}<$a[2]?$a[2]:$hxb{$a[1]}}
if(!exists $hya{$a[1]}){$hya{$a[1]}=$a[3];$hyb{$a[1]}=$a[3]}else{$hya{$a[1]}=$hya{$a[1]}>$a[3]?$a[3]:$hya{$a[1]};$hyb{$a[1]}=$hyb{$a[1]}<$a[3]?$a[3]:$hyb{$a[1]}}
}

for $i(sort {$a cmp $b} keys %hxa){
$xa=$hxa{$i}-int(($hxb{$i}-$hxa{$i})/2);
$xb=$hxb{$i}+int(($hxb{$i}-$hxa{$i})/2);
$xx=$xa+int(($xb-$xa)/2);
$ya=$hya{$i}-int(($hyb{$i}-$hya{$i})/2);
$yb=$hyb{$i}+int(($hyb{$i}-$hya{$i})/2);
$yy=$ya+int(($yb-$ya)/2);

print "$i\t$hxa{$i}\t$hxb{$i}\t$hya{$i}\t$hyb{$i}\t";

if(($xb-$xa)>($yb-$ya)){$ya=$yy-int(($xb-$xa)/2);$yb=$yy+int(($xb-$xa)/2);}
if(($xb-$xa)<($yb-$ya)){$xa=$xx-int(($yb-$ya)/2);$xb=$xx+int(($yb-$ya)/2);}

print "$xa\t$xb\t$ya\t$yb\n";

}
