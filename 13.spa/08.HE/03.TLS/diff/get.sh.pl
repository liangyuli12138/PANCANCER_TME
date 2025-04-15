while(<>){
chomp;@a=split;
$x=$a[1].".".$a[0];
$o=`cat findmarker.py`;
$o=~s/aaaa/$a[0]/g;
open OUT0,">shell/$x.findmarker.py";print OUT0 $o;

$z="/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/08.HE/03.TLS/at/$a[0].at";
undef %hb;
open IN,"$z";<IN>;while(<IN>){chomp;@b=split(/,/);$hb{$b[1]}=1};
for $i(keys %hb){
if($i eq "others"){next};
#print "shell/$x.$i.sh\n";
$y="shell/$x.$i.sh";
open OUT,">$y";
print OUT "/jdfssz1/ST_TSCBI/P22Z10200N0433/USER/zhangzhao/software/anaconda3/bin/python /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/08.HE/03.TLS/diff/shell/$x.findmarker.py /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/$a[0]","_cellbin.final.celltype.h5ad $i /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/08.HE/03.TLS/diff/out/$x/$x.$i.diff.csv\n"
}}
