open IN0,$ARGV[0];
<IN0>;
while(<IN0>){
chomp;s/\"//g;
@a=split(/,/);
$ha{$a[7]}=$a[-3];$haa{$a[1]}=1;$haaa{$a[2]}=1;

$hbc{$a[1]}{$a[7]}++;
$hb{$a[1]}{$a[7]}++;
$hbb{$a[1]}{$a[7]}+=$a[5];

$hde{$a[2]}{$a[7]}++;
$hd{$a[2]}{$a[7]}++;
$hdd{$a[2]}{$a[7]}+=$a[5];
};

open IN1,$ARGV[1];
<IN1>;
while(<IN1>){
chomp;s/\"//g;
@a=split(/,/);
$ha{$a[7]}=$a[-3];$haa{$a[1]}=1;$haaa{$a[2]}=1;

$hbc{$a[1]}{$a[7]}++;
$hc{$a[1]}{$a[7]}++;
$hcc{$a[1]}{$a[7]}+=$a[5];

$hde{$a[2]}{$a[7]}++;
$he{$a[2]}{$a[7]}++;
$hee{$a[2]}{$a[7]}+=$a[5];
};

open IN2,$ARGV[2];
<IN2>;
while(<IN2>){
chomp;s/\"//g;
@a=split(/,/);
$hg{$a[0]}="$a[1]\t$a[2]";
}

print "source_target\tcelltype\tinteraction_name\tCCI-Lymphoid_1_2_num\tCCI-Lymphoid3_num\ttotal\tadvantage_diff%\tadvantage\tCCI-Lymphoid_1_2_prob\tCCI-Lymphoid3_prob\tpathway_name\tgene1\tbaseMean1\tlog2FoldChange1\tgene2\tbaseMean2\tlog2FoldChange2\tgene3\tbaseMean3\tlog2FoldChange3\n";

for $i(sort {$a cmp $b} keys %ha){
for $j(sort {$a cmp $b} keys %haa){
$hb{$j}{$i}=$hb{$j}{$i}?$hb{$j}{$i}:0;
$hc{$j}{$i}=$hc{$j}{$i}?$hc{$j}{$i}:0;
$hbb{$j}{$i}=$hbb{$j}{$i}?$hbb{$j}{$i}:0;
$hcc{$j}{$i}=$hcc{$j}{$i}?$hcc{$j}{$i}:0;
if($hb{$j}{$i}==0 && $hc{$j}{$i}==0){next};
if($hbc{$j}{$i}==0){$pbc=0}else{$pbc=(abs($hb{$j}{$i}-$hc{$j}{$i}))/$hbc{$j}{$i}};
$o=$hc{$j}{$i}>$hb{$j}{$i}?"Lymphoid3":"Lymphoid_1_2";
@a=split(/\_/,$i);
$ho{$p} .= "source\t$j\t$i\t$hb{$j}{$i}\t$hc{$j}{$i}\t$hbc{$j}{$i}\t$pbc\t$o\t$hbb{$j}{$i}\t$hcc{$j}{$i}\t$ha{$i}";
for($j=0;$j<@a;$j++){$hg{$a[$j]}=$hg{$a[$j]}?$hg{$a[$j]}:"NA\tNA";$ho{$p} .=  "\t$a[$j]\t$hg{$a[$j]}"};
$ho{$p} .=  "\n"

}
for $j(sort {$a cmp $b} keys %haaa){
$hd{$j}{$i}=$hd{$j}{$i}?$hd{$j}{$i}:0;
$he{$j}{$i}=$he{$j}{$i}?$he{$j}{$i}:0;
$hdd{$j}{$i}=$hdd{$j}{$i}?$hdd{$j}{$i}:0;
$hee{$j}{$i}=$hee{$j}{$i}?$hee{$j}{$i}:0;
if($hd{$j}{$i}==0 && $he{$j}{$i}==0){next};
if($hde{$j}{$i}==0){$pde=0}else{$pde=(abs($hd{$j}{$i}-$he{$j}{$i}))/$hde{$j}{$i}};
$o=$he{$j}{$i}>$hd{$j}{$i}?"Lymphoid3":"Lymphoid_1_2";
@a=split(/\_/,$i);
$ho{$p} .= "target\t$j\t$i\t$hd{$j}{$i}\t$he{$j}{$i}\t$hde{$j}{$i}\t$pde\t$o\t$hdd{$j}{$i}\t$hee{$j}{$i}\t$ha{$i}";
for($j=0;$j<@a;$j++){$hg{$a[$j]}=$hg{$a[$j]}?$hg{$a[$j]}:"NA\tNA";$ho{$p} .=  "\t$a[$j]\t$hg{$a[$j]}"};
$ho{$p} .=  "\n"

}}

for $i(sort {$b <=> $a} keys %ho){print "$ho{$i}"}

