open IN0,$ARGV[0];<IN0>;
while(<IN0>){chomp;@a=split;if($a[1] eq "other"){next};$ha{$a[0]}=$a[1]};

open IN1,$ARGV[1];
while(<IN1>){chomp;$hb{$_}=1}

open IN3,$ARGV[3];
while(<IN3>){chomp;$hd{$_}=1}

open IN2,$ARGV[2];
while(<IN2>){
chomp;
@a=split(/,/);
$b=$a[-3];$b=~s/Cluster\d+//;$b=$b."_".$a[-2];
if(exists $hb{$b} && exists $ha{$a[-7]} && exists $hd{$a[-3]}){
$loc=$1 if($a[0]=~/(\d+\_\d+$)/);
$c=$a[-3]."_".$a[-2];
$ho{$c} .= "$loc,$ha{$a[-7]},$a[-7],$a[-3],$a[-2],$a[-1]\n"
}
}

for $i(keys %ho){
open IN,">at/$i.at";
print IN "$ho{$i}";
}
