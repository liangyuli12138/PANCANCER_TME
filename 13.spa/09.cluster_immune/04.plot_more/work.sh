cut -f 1 sn.list |while read i;do echo perl get.at.pl ../03.net_more/out100/$i.output.txt /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/$i""_cellbin.final.celltype.obs.csv \> at100/$i.immune.at;done|sh &

cut -f 1 sn.list |while read i;do echo perl get.at2.pl ../03.net_more/out100/$i.output.ex.txt /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/$i""_cellbin.final.celltype.obs.csv $i;done|sh &

cut -f 1 sn.list |while read i;do echo perl check.pl ../04.plot/at/$i.ori.at at100/$i.ori.at;done|sh |les

