open IN0,$ARGV[0];
<IN0>;
$t=$ARGV[1];
while(<IN0>){
chomp;
@a=split(/,/);
if($a[5]=~/(Lymphoid)/){$l="Lymphoid"}else{if($a[5]=~/(Myeloid)/){$l="Myeloid"}};
if($a[5]=~/Lymphoid/ || $a[5]=~/Myeloid/){
$ha{$t}{$l}++;
#print "$t\t$l\n";
if($a[3]=~/Cluster/){$ha{$l}{"aggregation"}++;if($a[4] =~ /malignant/){$ha{"aggregation"}{"malignant"}++}else{$ha{"aggregation"}{"normal"}++}}
else{$ha{$l}{"dispersion"}++;if($a[4] =~ /malignant/){$ha{"dispersion"}{"malignant"}++}else{$ha{"dispersion"}{"normal"}++}}

}
}

for $i(keys %ha){for $j(keys %{$ha{$i}}){print "$t.$i,$t.$j,$ha{$i}{$j}\n"}}
