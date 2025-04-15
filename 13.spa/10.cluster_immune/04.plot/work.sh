cut -f 1 sn.list |while read i;do echo perl get.at.pl immune.list  /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/$i""_cellbin.final.celltype.obs.csv /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/08.HE/01.Hexcell/region_merge/$i""_region.csv \> at/$i.all.at;done|sh &

cut -f 1 all.sn.list |while read i;do echo perl get.at.pl immune.list  /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/$i""_cellbin.final.celltype.obs.csv /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/08.HE/01.Hexcell/region_merge/$i""_region.csv \> at/$i.all.at;done|sh &

cut -f 1 all.sn.list |while read i;do echo perl get.at.icar.pl immune.list /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/09.cluster_immune/04.plot/at/$i"".ori.at /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/08.HE/01.Hexcell/region_merge/$i""_region.csv at/$i"".all.at \> at_icar/$i.icar.at;done|sh &

perl -e 'while(<>){chomp;@a=split;print "perl get.data.pl at_icar/$a[0].icar.at $a[1]\n"}' sn.list|sh > sample.for.sankey.csv

perl -e 'while(<>){chomp;@a=split(/\,/);$ha{$a[0]}{$a[1]}+=$a[2]};for $i(sort {$a cmp $b} keys %ha){for $j(sort {$a cmp $b} keys %{$ha{$i}}){print "$i,$j,$ha{$i}{$j}\n"}}'  sample.for.sankey.csv > tissue.for.sankey.csv

perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;@a=split(/,/);$ha{$a[0]}=$a[1]};open IN1,$ARGV[1];while(<IN1>){chomp;@a=split(/,/);$n=0;for $i(keys %ha){if($a[1]=~/$i/){print "$a[1],$ha{$i}\n"}}};open IN2,$ARGV[2];while(<IN2>){chomp;@a=split;print "$a[1].$a[1],#AEC7E8\n"}' color.csv tissue.for.sankey.csv sn.list > tissue.color.csv


