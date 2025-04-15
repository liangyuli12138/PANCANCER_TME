$s=$ARGV[3];chomp $s;

open IN0,$ARGV[0];
<IN0>;
while(<IN0>){
#f(!/$s/){next};
chomp;s/\"//g;
@a=split(/,/);
if($a[2]!~/$s/){next};
$ha{$a[7]}++;
$hb{$a[7]}++;};

open IN1,$ARGV[1];
<IN1>;
while(<IN1>){
#if(!/$s/){next};
chomp;s/\"//g;
@a=split(/,/);
if($a[2]!~/$s/){next};
$ha{$a[7]}++;
$hc{$a[7]}++};

open IN2,$ARGV[2];
<IN2>;
while(<IN2>){
chomp;s/\"//g;
@a=split(/,/);
$hg{$a[0]}="$a[1]\t$a[2]";
}

print "interaction_name\tCCI-Lymphoid_1_2\tCCI-Lymphoid3\ttotal\tadvantage_diff%\tadvantage\tgene1\tbaseMean1\tlog2FoldChange1\tgene2\tbaseMean2\tlog2FoldChange2\tgene3\tbaseMean3\tlog2FoldChange3\n";

for $i(sort {$a cmp $b} keys %ha){
$p=(abs($hc{$i}-$hb{$i}))/$ha{$i};$o=$hc{$i}>$hb{$i}?"Lymphoid3":"Lymphoid_1_2";
@a=split(/\_/,$i);$hb{$i}=$hb{$i}?$hb{$i}:0;$hc{$i}=$hc{$i}?$hc{$i}:0;
$ho{$p} .= "$i\t$hb{$i}\t$hc{$i}\t$ha{$i}\t$p\t$o";
for($j=0;$j<@a;$j++){$hg{$a[$j]}=$hg{$a[$j]}?$hg{$a[$j]}:"NA\tNA";$ho{$p} .=  "\t$a[$j]\t$hg{$a[$j]}"};
$ho{$p} .=  "\n"
}

for $i(sort {$b <=> $a} keys %ho){print "$ho{$i}"}

