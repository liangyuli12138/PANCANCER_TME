open IN0,$ARGV[0];
while(<IN0>){
chomp;
$num=$1 if(/Size\:\s+(\d+)/);
$id=$1 if(/Cluster\s+(\d+)/);
$id="Cluster"."$id";
s/^.+\)\:\s+//;
s/\s+//g;
@a=split(/\,/);
undef %hx;undef %hy;
for($i=0;$i<@a;$i++){
@b=split(/\_/,$a[$i]);
$hx{$b[0]}=1;
$hy{$b[1]}=1;
#print "$b[0]\t$b[1]\n";
$ha{$a[$i]}=$id;
}
for $i(sort {$a <=> $b} keys %hx){if(!exists $hnx{$id}){$hnx{$id}=$i-200};$hmx{$id}=$i+200}
for $i(sort {$a <=> $b} keys %hy){if(!exists $hny{$id}){$hny{$id}=$i-200};$hmy{$id}=$i+200}
}

print "cell,Cluster_net,Cluster_cell\n";
open IN1,$ARGV[1];
<IN1>;
while(<IN1>){
chomp;
chomp;@b=split(/,/);
if(exists $ha{$b[0]}){print "$b[0],$ha{$b[0]}"}else{print "$b[0],other"};
$cc="other";
for $i(keys %hnx){
if($b[-5]>=$hnx{$i} && $b[-5]<=$hmx{$i} && $b[-4]>=$hny{$i} && $b[-4]<=$hmy{$i}){$cc=$b[-3]}
}
print ",$cc\n"
}
