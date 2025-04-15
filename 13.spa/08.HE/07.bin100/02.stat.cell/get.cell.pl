open IN0,$ARGV[0];
<IN0>;
while(<IN0>){
chomp;
@a=split(/,/);$ha{$a[0]}=$a[1]};

open IN1,$ARGV[1];
$z=<IN1>;chomp $z;
@t=split(/,/,$z);

print "sample,variable,score\n";
while(<IN1>){
chomp;
@a=split(/,/);
if(!exists $ha{$a[0]}){next};
for($i=1;$i<@a;$i++){print "$t[$i],Cluster$ha{$a[0]},$a[$i]\n"}
}

