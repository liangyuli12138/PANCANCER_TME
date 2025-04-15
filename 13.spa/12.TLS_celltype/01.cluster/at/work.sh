perl -e '<>;while(<>){chomp;@a=split(/,/);if($a[0]=~/(\S+)\_(\d+\_\d+)$/){$x=$1;$y=$2};$ha{$x}{$y}="$a[1],$a[2],$a[3]";};for $i(keys %ha){$in="at1/$i.tls.at";open IN,$in;<IN>;$o="at/$i.tls.at";open OUT,">$o";print OUT "cell,celltype,cluster,new_groups,LM,allcluster\n";while(<IN>){chomp;@a=split(/,/);print OUT "$_,$ha{$i}{$a[0]}\n"}}' merge_25chip_immune_area.obs.new.at

