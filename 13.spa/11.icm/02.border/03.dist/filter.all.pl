open IN0,$ARGV[0];
while(<IN0>){chomp;@a=split;$ha{$a[0]}=$a[1]}

open IN1,$ARGV[1];
while(<IN1>){
chomp;
@x=split(/\//);$c="out/$x[-1].out";
$a="allcell/$x[-1]";$a=~s/all\.at\S+/all\.at.all/;
$b=$x[-1];$b=~s/\.all\.at\S+//;
open IN,$c;open OUT,">>$a";
while(<IN>){chomp;@f=split(/\t/);if($f[1]=~/dddd/){}else{
if($f[-1]>=$f[-2]){$o=$f[-1]*0.5;$o=$o<15?0:$o-15;}else{$o=-$f[-2]*0.5;$o=$o>-15?0:$o+15};
print OUT "$ha{$b},$b,$f[0],$f[1],$o\n";
}}
}
