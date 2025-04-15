cut -d "," -f 1,5  all_cell_meta.csv > all_cell_meta.csv.at

perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;@a=split(/,/);if($a[1] eq "Lymphoid1" || $a[1] eq "Lymphoid2"){$s="early"}elsif($a[1] eq "Lymphoid3"){$s="late"}else{$s="other"};$ha{$a[0]}="$a[1],$s"};open IN1,$ARGV[1];while(<IN1>){chomp;@a=split(/,/);$x=$a[1].$a[3];$y=$a[1]."_".$a[2];if(exists $ha{$x}){$hb{$y}="$ha{$x},$a[3],$a[4]"}};print "cell,new_groups,stage,cluster,celltype\n";open IN2,$ARGV[2];<IN2>;while(<IN2>){chomp;@a=split(/,/);if(exists $hb{$a[0]}){print "$a[0],$hb{$a[0]}\n"}}' l2m.stat.merge.csv all_cell_meta.csv merge_25chip_immune_area.obs > all_cell_meta.csv.at

