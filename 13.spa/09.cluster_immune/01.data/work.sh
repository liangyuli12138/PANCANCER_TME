for i in /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/*cellbin.final.celltype.obs.csv;do echo perl get.at.pl $i immune.list ;done|sh

perl -e 'while(<>){chomp;$nu=$1 if(/Size: (\d+)/);if($nu>=100){print "$_\n"}}' output.txt > output.txt.filter

