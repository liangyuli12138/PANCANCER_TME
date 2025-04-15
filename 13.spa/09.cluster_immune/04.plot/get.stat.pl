$t=$ARGV[0];
$ii=$1 if($t=~/out\/(\S+)\.output/);
open IN0,$ARGV[0];
while(<IN0>){
chomp;if(/Boundary/){next};
$num=$1 if(/Size\:\s+(\d+)/);
$id=$1 if(/Cluster\s+(\d+)/);
$total+=$num;
$xx=$id
}

print "$ii\t$xx\t$total\n";
