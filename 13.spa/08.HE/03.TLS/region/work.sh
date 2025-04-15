perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;$t=$_;$t=~s/_immue_region.csv//;open IN,$_;<IN>;while(<IN>){chomp;@a=split(/,/);if(!exists $hc{$t}{$a[1]}){$hb{$t}++};$hc{$t}{$a[1]}=1;$ha{$t}++}};open IN1,$ARGV[1];while(<IN1>){chomp;@a=split;$ha{$a[0]}=$ha{$a[0]}?$ha{$a[0]}:0;$hb{$a[0]}=$hb{$a[0]}?$hb{$a[0]}:0;$y+=$hb{$a[0]};$z+=$ha{$a[0]};print "$a[1]\t$a[0]\t$hb{$a[0]}\t$ha{$a[0]}\n"};print "total\t$y\t$z\n"' region.list sn.list > tls.region.list.stat.xls

cat *immue_region.csv|perl get.window.pl > all.tls.windows.list

perl -e 'while(<>){chomp;@a=split;print "perl stat.pl merge.list /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/08.HE/03.TLS/region/$a[0]_immue_region.csv /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/$a[0]_cellbin.final.celltype.obs.csv\n"}' sn.list |sh > merge.stat

perl -e 'while(<>){chomp;if(/Unk/){next};@a=split;$hb{$a[1]}=1;$ha{$a[0]}{$a[1]}=$a[2]};for $i(sort {$a cmp $b} keys %hb){print "\t$i"};print "\n";for $j(sort {$a cmp $b} keys %ha){print "$j";for $i(sort {$a cmp $b} keys %hb){$ha{$j}{$i}=$ha{$j}{$i}?$ha{$j}{$i}:0;print "\t$ha{$j}{$i}"};print "\n"}' merge.stat > merge.stat.xls

perl -e 'while(<>){chomp;@a=split;print "perl stat.pl merge.list /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/08.HE/03.TLS/region/$a[0]_immue_region.csv /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/$a[0]_cellbin.final.celltype.obs.csv\n"}' sn.list |sh > merge.stat.sub

