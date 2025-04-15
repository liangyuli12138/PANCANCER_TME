while(<>)
{chomp;$t=$_;$t=~s/\.out\S+\.csv//;
open IN,$_;
$p=<IN>;chomp $p;@a=split(/,/,$p);for($i=1;$i<@a;$i++){$ha{$a[$i]}=1};

while(<IN>)
{chomp;@b=split(/,/);
for($j=1;$j<@b;$j++){
$hb{$t}{$b[0]}{$a[$j]}=$b[$j]
}
}
}

for $k(sort {$a cmp $b} keys %hb){
print "tissue\tgruops";
for $i(sort {$a cmp $b} keys %ha){print "\t$i"};print "\n";
for $j(sort {$a cmp $b} keys %{$hb{$k}}){print "$k\t$j";
for $i(sort {$a cmp $b} keys %ha){print "\t$hb{$k}{$j}{$i}"};print "\n";
}
print "\n"
}
