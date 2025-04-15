open IN2,$ARGV[2];
while(<IN2>){
chomp;
@a=split(/,/);
$hm{$a[0]}=1
}

open IN0,$ARGV[0];
while(<IN0>){
chomp;
$ha{$_}=1}

print "\"\",\"id\",\"celltype\"\n";
open IN1,$ARGV[1];
<IN1>;
while(<IN1>){
chomp;$n++;
chomp;@b=split(/,/);
if(exists $hm{$b[0]}){}else{print "\"$n\",\"$b[0]\",\"check\"\n"};

}
