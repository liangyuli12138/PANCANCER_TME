open IN0,$ARGV[0];$t=$ARGV[0];<IN0>;while(<IN0>){chomp;$n++;if($n<=100){@a=split(/,/);$ha{$a[1]}=$a[0];$haa{$a[1]}=$_}};
open IN1,$ARGV[1];$p=$ARGV[1];<IN1>;while(<IN1>){chomp;$m++;if($m<=100){@a=split(/,/);$hb{$a[1]}=$a[0];$hbb{$a[1]}=$_}};
for $i(keys %ha){if(exists $hb{$i}){$x=$ha{$i}+$hb{$i};$hc{$x}{$i}=$haa{$i};$hcc{$x}{$i}=$hbb{$i};}}
for $i(sort {$a <=> $b} keys %hc){for $j(keys %{$hc{$i}}){$y++;print "$y,$j,$hc{$i}{$j}\n$y,$j,$hcc{$i}{$j}\n"}}
