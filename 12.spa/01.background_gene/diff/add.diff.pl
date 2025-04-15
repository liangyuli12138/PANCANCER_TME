open IN0,$ARGV[0];
while(<IN0>){
chomp;@a=split(/\//);$a[-1]=~s/diff\.//;$a[-1]=~s/\.csv//;$n=0;open IN,$_;while(<IN>){chomp;@b=split(/\,/);$n++;if($n>100){last};$b[-2]=sprintf("%.2f",$b[-2]);$b[-1]=sprintf("%.2f",$b[-1]);$ha{$b[1]}.="{$a[-1]|$b[-2]|$b[-1]};"}};
open IN1,$ARGV[1];
while(<IN1>){chomp;
@a=split;
$ha{$a[0]}=$ha{$a[0]}?$ha{$a[0]}:"NA";
print "$_\t$ha{$a[0]}\n";
}
