open IN3,$ARGV[3];
<IN3>;
while(<IN3>){
chomp;
s/\"//g;
@a=split(/,/);
$hc{$a[1]}=$a[2]
}

open IN2,$ARGV[2];
while(<IN2>){
chomp;
@a=split(/,/);
if($a[1] eq "malignant"){$hm{$a[0]}="malignant"};
}

open IN0,$ARGV[0];
while(<IN0>){
chomp;
$ha{$_}=1}

print "\"\",\"id\",\"celltype\",\"cluster\",\"malignant\",\"immune\"\n";
open IN1,$ARGV[1];
<IN1>;
while(<IN1>){
s/\"//g;
chomp;$n++;
chomp;@b=split(/,/);
if(exists $ha{$b[2]}){$m=$hm{$b[1]}?$hm{$b[1]}:"other";print "\"$n\",\"$b[1]\",\"$b[2]\",\"$b[3]\",\"$m\",\"$hc{$b[1]}\"\n"}elsif(exists $hm{$b[1]}){print "\"$n\",\"$b[1]\",\"malignant\",\"malignant\",\"malignant\",\"$hc{$b[1]}\"\n"}else{print "\"$n\",\"$b[1]\",\"other\",\"other\",\"other\",\"$hc{$b[1]}\"\n"};

}
