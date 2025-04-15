open IN0,$ARGV[0];
<IN0>;
while(<IN0>){
chomp;s/\"//g;
@a=split(/,/);
if($a[2]=~/malignant/){$ha{$a[1]}=1}
}

open IN1,$ARGV[1];
$t=<IN1>;print "$t";
while(<IN1>){
chomp;
@a=split(/,/);
$a[1]=~s/\"//g;
if(exists $ha{$a[1]}){s/other/malignant/g}
print "$_\n"
}
