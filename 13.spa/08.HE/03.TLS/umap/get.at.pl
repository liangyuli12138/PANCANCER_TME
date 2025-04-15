open IN0,$ARGV[0];
while(<IN0>){chomp;@a=split;$ha{$a[0]}=$a[1]};

open IN1,$ARGV[1];
while(<IN1>){chomp;@a=split(/,/);$hb{$a[1]}{$a[0]}=1;$a[1]=~s/\_TLS\S+//;$hc{$a[1]}{$a[0]}=1;}

open IN2,$ARGV[2];
while(<IN2>){
chomp;$t=$_;
$id=$1 if(/data\/(\S+)\_cellbin/);
open IN,$t;
while(<IN>){chomp;@a=split(/,/);$hd{$id}{$a[0]}=$ha{$a[-3]}}
}

open IN3,$ARGV[3];
while(<IN3>){
chomp;
@a=split;$o="$a[0].at";
open OUT,">at/$o";
print OUT ",id,celltype\n";
$t=$a[0];$t=~s/\_TLS\S+//;
for $i(keys %{$hd{$t}}){
@b=split(/\_/,$i);
if($b[0]>$a[5] && $b[0]<$a[6] && $b[1]>$a[7] && $b[1]<$a[8]){
$n++;
$o=$hb{$a[0]}{$i}?$hd{$t}{$i}:"others";
print OUT "$n,$i,$o\n"
}
}}
