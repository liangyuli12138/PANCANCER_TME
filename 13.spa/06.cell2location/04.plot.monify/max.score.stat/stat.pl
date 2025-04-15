@t=(0,0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1);
$m="m1";
while(<>){
chomp;
if(/max_score/){next};
@a=split(/,/);$hb{$a[-3]}=1;$c=0;
for($i=0;$i<@t-1;$i++){$j=$i+1;if($a[-2]>$t[$i] && $a[-2]<=$t[$j]){$ha{$t[$j]}{$a[-3]}++;$c=1;last}};if($c==0){$ha{$m}{$a[-3]}++}};

print "index";
for $i(sort {$a cmp $b} keys %hb){print "\t$i"};
print "\n";
for ($j=1;$j<@t;$j++){print "$t[$j]";for $i(sort {$a cmp $b} keys %hb){$ha{$t[$j]}{$i}=$ha{$t[$j]}{$i}?$ha{$t[$j]}{$i}:0;print "\t$ha{$t[$j]}{$i}"};print "\n"};
print "m1";
for $i(sort {$a cmp $b} keys %hb){$ha{$m}{$i}=$ha{$m}{$i}?$ha{$m}{$i}:0;print "\t$ha{$m}{$i}"};print "\n"
