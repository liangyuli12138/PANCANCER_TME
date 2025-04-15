open IN0,$ARGV[0];
while(<IN0>){
chomp;
@a=split;
$ha{$a[0]}=$a[2]}

print "id,region,x,y,celltype\n";
open IN1,$ARGV[1];
<IN1>;
while(<IN1>){
chomp;
for $i(keys %ha){s/$i/$ha{$i}/;}
print "$_\n"
}
