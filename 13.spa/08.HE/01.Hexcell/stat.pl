open IN2,$ARGV[2];while(<IN2>){chomp;@a=split;$hx{$a[0]}=$a[1]}

open IN0,$ARGV[0];
while(<IN0>){
chomp;
@a=split(/\//);
@b=split(/\./,$a[-1]);
$b[1]=~s/\_meta\_ori//;
$b[1]=~s/\_meta\_ori//;
open IN,$_;<IN>;
while(<IN>){
chomp;s/\"//g;@c=split(/,/);
@d=split(/\_/,$c[2]);
$ha{$b[1]}{$c[1]}=$c[2];
$hb{$b[1]}{$c[1]}=$d[0];
$haa{$c[2]}=1;$hbb{$d[0]}=1;
}};

open IN1,$ARGV[1];
while(<IN1>){
chomp;
@a=split(/\//);
@b=split(/\_/,$a[-1]);
open IN,$_;<IN>;
while(<IN>){
chomp;
@c=split(/,/);
$hc{$b[0]}{$c[1]}{$ha{$b[0]}{$c[0]}}++;
$hd{$b[0]}{$c[1]}{$hb{$b[0]}{$c[0]}}++;
$hcc{$c[1]}=1;
}}

for $k(sort {$hx{$b} cmp $hx{$a}} keys %hx){

print "tissue\tid\tregion";
for $i(sort {$a cmp $b} keys %hbb){
print "\t$i"
};print "\n";

for $j(sort {$a cmp $b} keys %{$hd{$k}}){
print "$hx{$k}\t$k\t$j";
for $i(sort {$a cmp $b} keys %hbb){
$hd{$k}{$j}{$i}=$hd{$k}{$j}{$i}?$hd{$k}{$j}{$i}:0;
print "\t$hd{$k}{$j}{$i}"};print "\n"};print "\n"}

for $k(sort {$hx{$b} cmp $hx{$a}} keys %hx){

print "tissue\tid\tregion";
for $i(sort {$a cmp $b} keys %haa){
print "\t$i"
};print "\n";

for $j(sort {$a cmp $b} keys %{$hc{$k}}){
print "$hx{$k}\t$k\t$j";
for $i(sort {$a cmp $b} keys %haa){
$hc{$k}{$j}{$i}=$hc{$k}{$j}{$i}?$hc{$k}{$j}{$i}:0;
print "\t$hc{$k}{$j}{$i}"};print "\n"};print "\n"}











