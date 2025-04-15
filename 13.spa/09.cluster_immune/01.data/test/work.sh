perl get.at.pl output.txt.filter /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/D01972D1_cellbin.final.celltype.obs.csv > output.at

perl -e 'print "\"\",\"id\",\"celltype\"\n";<>;while(<>){chomp;@a=split(/\,/);$n++;print "\"$n\",\"$a[0]\",\"$a[1]\"\n"}' output.at > output.at.net
perl -e 'print "\"\",\"id\",\"celltype\"\n";<>;while(<>){chomp;@a=split(/\,/);$n++;print "\"$n\",\"$a[0]\",\"$a[2]\"\n"}' output.at > output.at.cell

