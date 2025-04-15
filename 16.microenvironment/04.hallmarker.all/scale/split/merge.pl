open IN0,$ARGV[0];
while(<IN0>){chomp;$ha{$_}=1};

print "sample,variable,score\n";
open IN1,$ARGV[1];
while(<IN1>){
chomp;$t=$_;
$t=~s/\.UCell\.score\.csv//;
open IN,$_;
$x=<IN>;chomp $x;$x=~s/\"//g;$x=~s/\_UCell//g;
@y=split(/,/,$x);
while(<IN>){
chomp;@a=split(/,/);
for($i=1;$i<@a;$i++){if(exists $ha{$y[$i]}){print "$t,$y[$i],$a[$i]\n"}}
}}
