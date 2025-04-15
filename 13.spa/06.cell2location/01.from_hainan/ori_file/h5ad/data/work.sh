perl -e 'while(<>){chomp;@a=split;$x=$a[0]."_cellbin.final.celltype.obs.csv";$y=$a[1]."_cellbin.final.celltype.obs.csv";print "mv Stereo_seq/$x Stereo_seq/$y\n"}' name.list |sh

