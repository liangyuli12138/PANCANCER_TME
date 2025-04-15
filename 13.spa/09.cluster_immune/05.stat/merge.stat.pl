open IN0,$ARGV[0];
while(<IN0>){
chomp;
@a=split;$ha{$a[0]}=1
}

open IN1,$ARGV[1];
while(<IN1>){
chomp;
$t=$_;$id=$1 if($t=~/score\/(\S+)\.at/);
open IN,$t;
<IN>;
while(<IN>){
chomp;
@a=split(/,/);
$hb{$id}{$a[0]}="$a[2]\t$a[3]\t$a[4]\t$a[5]\t$a[6]";
}
}

#open IN3,$ARGV[3];


print "sample\tcell\tcluster\ttype\timmune\tantiTumor\tproTumor\tangiogenesis_fibrosis\tmalignant_properties\thallmarkerEMT\n";
open IN2,$ARGV[2];
while(<IN2>){
chomp;
$t=$_;$id=$1 if($t=~/at\/(\S+)\.ori/);
open IN,$t;
<IN>;
while(<IN>){
chomp;s/\"//g;@a=split(/,/);if($a[-1]=~/other/){next};
$o=$ha{$a[2]}?"Immune":"others";
print "$id\t$a[1]\t$a[-1]\t$a[2]\t$o\t$hb{$id}{$a[1]}\n";
}}
