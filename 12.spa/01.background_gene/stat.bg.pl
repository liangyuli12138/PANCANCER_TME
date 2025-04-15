while(<>){
if(/^#/){next};
@a=split;
if(!exists $hb{$a[5]}){$b++;$hb{$a[5]}=1}
if(!exists $ha{$a[0]}{$a[5]}){$hg{$a[0]}++};$ha{$a[0]}{$a[5]}+=$a[3];
$hc{$a[0]}+=$a[3];
}

for $i(keys %ha){for $j(sort {$ha{$i}{$b} <=> $ha{$i}{$a}} keys %{$ha{$i}}){$hm{$i}++;if($hm{$i}>=($hg{$i}*0.3) && $hm{$i}<=($hg{$i}*0.3+1)){$hn{$i}=$ha{$i}{$j}}}}

for $i(sort {$hc{$b} <=> $hc{$a}} keys %hc){
$p=$hc{$i}/$hg{$i};$x=$hg{$i}/$b;
print "$i\t$hc{$i}\t$hg{$i}\t$p\t$x\t$hn{$i}\n"
}

