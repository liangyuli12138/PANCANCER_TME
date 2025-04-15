open IN1,$ARGV[1];
while(<IN1>){
chomp;
$t=$_;
$ha{$t}=1;
}

open IN0,$ARGV[0];
$t=$ARGV[0];
$id=$1 if($t=~/data\/(\S+)\_cellbin/);

open OUT0,">at/$id.input";
open OUT1,">at/$id.at";
print OUT0 "cell\n";
print OUT1 "cell,region\n";

while(<IN0>){
chomp;@a=split(/,/);
if(exists $ha{$a[-3]})
{
print OUT0 "$a[0]\n";print OUT1 "$a[0],Immune\n";
}
}
