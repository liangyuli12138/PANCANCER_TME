open IN0,$ARGV[0];
<IN0>;
while(<IN0>){
chomp;
@a=split(/,/);
if($a[0]=~/(\S+)\_\d+\_\d+/){$id=$1};
$a[0]=~s/$id\_//;$ha{$id}{$a[0]}=$a[2]};

open IN1,$ARGV[1];
while(<IN1>){
chomp;
$ip=$1 if(/\/data\/(\S+)\_cellbin/);
$aa="at/$ip.at";$bb="at/$ip.input";open OUT0,">$aa";open OUT1,">$bb";
print OUT0 "cell,tlsclu,batch\n";print OUT1 "cell\n";
open IN,$_;<IN>;
while(<IN>){chomp;@c=split(/,/);if(exists $ha{$ip}{$c[0]}){print OUT0 "$c[0],$ha{$ip}{$c[0]},$ip\n";print OUT1 "$c[0]\n"}}}
