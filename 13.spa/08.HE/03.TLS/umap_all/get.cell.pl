open IN0,$ARGV[0];
while(<IN0>){chomp;@a=split;$ha{$a[0]}=$a[1]};

open IN1,$ARGV[1];
while(<IN1>){
chomp;
@a=split(/,/);
$hb{$a[0]}=1;
}

open OUT0,">$ARGV[3]";
open OUT1,">$ARGV[4]";

print OUT0 ",id,celltype\n";
print OUT1 ",id,celltype\n";

open IN2,$ARGV[2];
<IN2>;
while(<IN2>){
#s/\"//g;
chomp;
@a=split(/,/);
if(exists $hb{$a[0]}){$n++;print OUT0 "$n,$a[0],$a[-3]\n";print OUT1 "$n,$a[0],$ha{$a[-3]}\n"}
else{$n++;print OUT0 "$n,$a[0],others\n";print OUT1 "$n,$a[0],others\n"}
}
