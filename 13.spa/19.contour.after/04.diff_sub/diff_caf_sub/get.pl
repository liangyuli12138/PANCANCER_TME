while(<>){
chomp;
@a=split;
$ha{$a[1]}=1};

@t=("0","0-100","100-300","300-500");

for($j=0;$j<@t;$j++){
for $i(keys %ha){
$x=$i."_Lymphoid_1_2_".$t[$j];
$y=$i."_Lymphoid3_".$t[$j];

print "/hwfssz4/BC_PUB/Software/07.User-defined/03.Animal_Plant/wubin/mamba/bin/python /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/19.contour.after/04.diff_sub/diff_caf/findmarker.py /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/19.contour.after/04.diff_sub/diff_caf/pancancer.icar.contour.cell.h5ad $x $y /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/19.contour.after/04.diff_sub/diff_caf/diff/diff.$x.csv\n/hwfssz4/BC_PUB/Software/07.User-defined/03.Animal_Plant/wubin/mamba/bin/python /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/19.contour.after/04.diff_sub/diff_caf/findmarker.py /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/19.contour.after/04.diff_sub/diff_caf/pancancer.icar.contour.cell.h5ad $y $x /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/19.contour.after/04.diff_sub/diff_caf/diff/diff.$y.csv\n"}}
