open IN2,$ARGV[2];
while(<IN2>){
chomp;
@a=split(/,/);
$hr{$a[0]}=1;
if($a[1] eq "malignant"){$hm{$a[0]}=1};
}

open IN0,$ARGV[0];
while(<IN0>){
chomp;
$ha{$_}=1}

print "\"\",\"id\",\"celltype\"\n";
open IN1,$ARGV[1];
<IN1>;
while(<IN1>){
chomp;@b=split(/,/);if(!exists $hr{$b[0]}){next};$n++;
if(exists $ha{$b[-3]}){print "\"$n\",\"$b[0]\",\"$b[-3]\"\n"}elsif(exists $hm{$b[0]}){print "\"$n\",\"$b[0]\",\"malignant\"\n"}else{print "\"$n\",\"$b[0]\",\"other\"\n"};

}
