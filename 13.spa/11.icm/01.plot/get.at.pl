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
if($a[2] eq "malignant" || $a[2] eq "other" ){print "$o\n"}
elsif(exists $ha{$a[2]}){$o=~s/$a[2]/$ha{$a[2]}/;print "$o\n"}
else{$o=~s/$a[2]/other/;print "$o\n"}
}
