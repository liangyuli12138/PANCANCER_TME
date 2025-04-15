perl -e 'print "cell\n";while(<>){chomp;@a=split(/,/);if($a[-1]=~/Lymphoid_B/ && $a[-2]=~/Lymphoid/){print "$a[0]\n"}}' pancancer.icar.all.cell.obs > filter.B.input
perl -e 'print "cell\n";while(<>){chomp;@a=split(/,/);if($a[-1]=~/Lymphoid_CD/ && $a[-2]=~/Lymphoid/){print "$a[0]\n"}}' pancancer.icar.all.cell.obs > filter.T.input
perl -e 'print "cell\n";while(<>){chomp;@a=split(/,/);if($a[-1]=~/Myeloid_cDC/ && $a[-2]=~/Lymphoid/){print "$a[0]\n"}}' pancancer.icar.all.cell.obs > filter.cDC.input

perl -e '@a=("pancancer.icar.B.cell.h5ad","pancancer.icar.cDC.cell.h5ad","pancancer.icar.T.cell.h5ad");@b=("Lymphoid_1_2","Lymphoid3");for($i=0;$i<@a;$i++){for($j=0;$j<@b;$j++){print "/hwfssz4/BC_PUB/Software/07.User-defined/03.Animal_Plant/wubin/mamba/bin/python /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/15.fig6/05.diff/B/findmarker.py /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/15.fig6/05.diff/B/$a[$i] $b[$j] /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/15.fig6/05.diff/B/diff/diff.$a[$i].$b[$j].csv\n"}}' > all.diff.sh

