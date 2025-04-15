perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;@a=split(/,/);$ha{$a[3]}=$a[1]};print ",sample,cell,cluster,type,new_groups,new_cluster\n";open IN1,$ARGV[1];<IN1>;while(<IN1>){chomp;@a=split(/,/);$t=$a[1].$a[3];if(exists $ha{$t}){print "$a[0],$a[1],$a[2],$a[3],$a[4],$ha{$t},$t\n"}}' merge_25chip_immune_area.obs.new.at.celltype all_cell_meta.csv > icar.all_cell_meta.csv

