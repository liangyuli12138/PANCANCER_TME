cut -f 1 all.sn.list |while read i;do echo perl get.at.pl immune.list  /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/$i""_cellbin.final.celltype.obs.csv /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/08.HE/01.Hexcell/region_merge/$i""_region.csv \> at.ori/$i.all.at;done|sh &
cut -f 1 all.sn.list |while read i;do echo perl get.icm.pl icm.list at.ori/$i.all.at \> at.icm/$i.all.at;done|sh

