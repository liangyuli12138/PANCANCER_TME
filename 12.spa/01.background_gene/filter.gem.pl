open IN0,$ARGV[0];
while(<IN0>){
chomp;
@a=split;
$ha{$a[0]}=$a[5]
}

print "geneID\tx\ty\tMIDCount\tExonCount\tlabel\ttag\n";
open IN1,$ARGV[1];
while(<IN1>){
if(/^#/ || /^geneID/){next};
chomp;
@a=split;if($a[5] eq "0"){next};
if(!exists $ha{$a[0]}){print "$a[0]\t$a[1]\t$a[2]\t$a[3]\t$a[4]\t$a[5]\t1\n"}else{
if(!exists $hb{$a[0]}{$a[5]}){$hb{$a[0]}{$a[5]}=$ha{$a[0]}}
if($a[3]>$hb{$a[0]}{$a[5]}){$a[3]=$a[3]-$hb{$a[0]}{$a[5]};$a[4]=$a[3];$hb{$a[0]}{$a[5]}=0;print "$a[0]\t$a[1]\t$a[2]\t$a[3]\t$a[4]\t$a[5]\t1\n"}else
{$hb{$a[0]}{$a[5]}=$hb{$a[0]}{$a[5]}-$a[3]}
}}
