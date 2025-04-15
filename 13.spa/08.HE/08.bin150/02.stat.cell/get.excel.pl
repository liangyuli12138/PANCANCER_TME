open IN0,$ARGV[0];
<IN0>;
while(<IN0>){
chomp;
@a=split(/,/);$ha{$a[0]}=$a[1]};

open IN1,$ARGV[1];
$z=<IN1>;chomp $z;
@t=split(/,/,$z);

while(<IN1>){
chomp;
@a=split(/,/);
if(!exists $ha{$a[0]}){next};
for($i=1;$i<@a;$i++){$x="$ha{$a[0]}";$hb{$x}{$t[$i]}+=$a[$i];$hc{$x}{$t[$i]}++}
}

$z=~s/,/\t/g;
print "$z\n";
for $i(sort {$a <=> $b} keys %hb){
print "Cluster$i";
for ($j=1;$j<@t;$j++){
$o=$hc{$i}{$t[$j]}>0?$hb{$i}{$t[$j]}/$hc{$i}{$t[$j]}:0;

print "\t$o"

}
print "\n"
}
