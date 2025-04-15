open IN1,$ARGV[1];
while(<IN1>){
chomp;$t=$_;$id=$1 if($t=~/get\.sf\/(\S+)_doublet/);
open IN,$t,<IN>;while(<IN>){
chomp;s/\"//g;@a=split(/\t/);$hb{$id}{$a[11]}++;$hi{$a[11]}=1;$hx{$id}{$a[0]}=1}};

open IN0,$ARGV[0];
while(<IN0>){
chomp;$t=$_;$id=$1 if($t=~/Merge\_(\S+)\_meta/);
@i=split(/\./,$id);
$hn{$i[1]}=$i[0];
open IN,$t,<IN>;while(<IN>){s/CD16\./CD16\+/;s/CD56\./CD56\+/;
chomp;s/\"//g;@a=split(/,/);
if(!exists $hx{$i[1]}{$a[1]}){next};
$ha{$i[1]}{$a[-1]}++}};

open IN2,$ARGV[2];
while(<IN2>){
chomp;$t=$_;$id=$1 if($t=~/Merge\_(\S+)\_meta/);
@i=split(/\./,$id);
open IN,$t,<IN>;while(<IN>){s/CD16\./CD16\+/;s/CD56\./CD56\+/;
chomp;s/\"//g;@a=split(/,/);
#if(!exists $hx{$i[1]}{$a[1]}){next};
$hc{$i[1]}{$a[-1]}++}};

open IN3,$ARGV[3];
while(<IN3>){
chomp;$t=$_;$id=$1 if($t=~/Merge\_(\S+)\_meta/);
@i=split(/\./,$id);
open IN,$t,<IN>;while(<IN>){s/CD16\./CD16\+/;s/CD56\./CD56\+/;
chomp;s/\"//g;@a=split(/,/);
#if(!exists $hx{$i[1]}{$a[1]}){next};
$hd{$i[1]}{$a[-3]}++};$hi{$a[-3]}=1};


print "sn\t0tumore\ttype";for $i(sort {$a cmp $b} keys %hi){print "\t$i"};print "\n";

for $j(sort {$a cmp $b} keys %ha){print "$j\t$hn{$j}\tcell2location1";
for $i(sort {$a cmp $b} keys %hi){$ha{$j}{$i}=$ha{$j}{$i}?$ha{$j}{$i}:0;print "\t$ha{$j}{$i}"};print "\n"};

for $j(sort {$a cmp $b} keys %hc){print "$j\t$hn{$j}\tcell2location2";
for $i(sort {$a cmp $b} keys %hi){$hc{$j}{$i}=$hc{$j}{$i}?$hc{$j}{$i}:0;print "\t$hc{$j}{$i}"};print "\n"};

for $j(sort {$a cmp $b} keys %hd){print "$j\t$hn{$j}\tcell2location3";
for $i(sort {$a cmp $b} keys %hi){$hd{$j}{$i}=$hd{$j}{$i}?$hd{$j}{$i}:0;print "\t$hd{$j}{$i}"};print "\n"};

for $j(sort {$a cmp $b} keys %hb){print "$j\t$hn{$j}\tRCTD";
for $i(sort {$a cmp $b} keys %hi){$hb{$j}{$i}=$hb{$j}{$i}?$hb{$j}{$i}:0;print "\t$hb{$j}{$i}"};print "\n"};
