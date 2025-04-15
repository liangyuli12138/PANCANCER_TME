open IN0,$ARGV[0];
while(<IN0>){
chomp;@a=split;$ha{$a[0]}=$a[1]}

open IN1,$ARGV[1];
$t=<IN1>;print "$t";
while(<IN1>){
chomp;
$o=$_;
s/\"//g;
@a=split(/,/);
if($a[4] eq "malignant" || $a[4] eq "other" ){print "$o\n"}
elsif(exists $ha{$a[4]}){$o=~s/$a[4]/$ha{$a[4]}/;print "$o\n"}
else{$o=~s/$a[4]/other/;print "$o\n"}
}
