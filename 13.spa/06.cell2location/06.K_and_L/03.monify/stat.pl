<>;
while(<>){
chomp;
@a=split(/,/);
$ha{$a[-4]}{$a[-1]}++;
$hb{$a[-1]}=1;
}

print "group";
for $j(sort {$a <=> $b} keys %hb){print "\t$j"};print "\n";
for $i(sort {$a cmp $b} keys %ha){
print "$i";
for $j(sort {$a <=> $b} keys %hb){print "\t$ha{$i}{$j}"};
print "\n";
}
