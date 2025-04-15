print "icc,gene,cell,priority,fc_log2,pct_lym3,pct_lym12\n";
open IN0,$ARGV[0];
while(<IN0>){
chomp;
$id=$1 if(/diff\.(\S+)\_Lymphoid/);
$l=$1 if(/\_Lymphoid(\S+)\.csv/);

open IN,$_;
;<IN>;
while(<IN>){
chomp;if($l==3){
@a=split(/,/);$ha{$a[1]}{$id}="$a[0],$a[1],$a[3],$a[-1]";$hc{$a[1]}{$id}=$a[0]}else{@a=split(/,/);$hb{$a[1]}{$id}="$a[-1]"}}};open IN1,$ARGV[1];while(<IN1>){chomp;@a=split(/\_/);for($i=0;$i<@a;$i++){for $j(sort {$hc{$a[$i]}{$a} <=> $hc{$a[$i]}{$b}} keys %{$hc{$a[$i]}}){print "$_,$a[$i],$j,$ha{$a[$i]}{$j},$hb{$a[$i]}{$j}\n"}}}
