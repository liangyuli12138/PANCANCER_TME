$t=$ARGV[2];
chomp $t;
open IN0,$ARGV[0];
while(<IN0>){
chomp;
if(/$t/){s/$t//;@a=split(/\,/);
$ha{$a[0]}=$a[1]}};

open IN1,$ARGV[1];
$o=<IN1>;
chomp $o;
$o=~s/\"\"\,\"id\"\,\"celltype\"/\"cell\"\,\"new_celltype\"/;
print "$o,new_groups\n";
while(<IN1>){
chomp;
s/\"\d+\"\,//;

@a=split(/,/);
$a[-1]=~s/\"//g;
#$o=$ha{$a[-1]}?$ha{$a[-1]}:"\"other\"";
if(!exists $ha{$a[-1]}){
if($a[-1]=~/Cluster/){print "$_,Cluster_other\n"}
else{print "$_,Cluster_out\n"}
}
else{print "$_,$ha{$a[-1]}\n"}}
