cut -f 1 ../sn.list|grep B01615B3|while read i;do echo perl get.at.pl ../immune.list zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/$i""_cellbin.final.celltype.obs.csv /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/08.HE/01.Hexcell/region_merge/$i""_region.csv \> $i.all.at;done|sh

cut -f 1 ../sn.list|grep B01615B3|while read i;do echo perl get.check.pl ../immune.list /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/$i""_cellbin.final.celltype.obs.csv /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/08.HE/01.Hexcell/region_merge/$i""_region.csv \> $i.all.at;done|sh

