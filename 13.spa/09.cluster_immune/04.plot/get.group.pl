open IN0,$ARGV[0];
while(<IN0>){
chomp;
@a=split(/,/);
if($a[1] eq "malignant"){$hm{$a[0]}=1};
}

open IN1,$ARGV[1];
<IN1>;
while(<IN1>){
chomp;
@a=split(/,/);
$ha{$a[0]}=$a[1];
}


open IN2,$ARGV[2];
$t=<IN2>;print "$t";
$ip=$1 if($ARGV[2]=~/at\/(\S+)\.immune/);
while(<IN2>){
chomp;s/\"//g;
@a=split(/,/);
if(exists $hm{$a[1]} && $a[2] eq "other"){print "\"$a[0]\",\"$a[1]\",\"malignant\"\n";next};
if($a[2] eq "other"){print "\"$a[0]\",\"$a[1]\",\"other\"\n";next};
if($a[2] ne "other"){$id=$ip."_".$a[2];print "\"$a[0]\",\"$a[1]\",\"group$ha{$id}\"\n"}
}
