$tti=$ARGV[0];
@a=split(/,/,$tti);
$tt=$a[0];$pp=$a[1];

open OUT0,">/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/05.diff_six/diff/find_diff/find_plot/$tt/find.$tt.py";

print OUT0 "import os
import datetime
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context

os.chdir(\"/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/05.diff_six/diff/find_diff/find_plot/$tt\");
adata_concat = sc.read_h5ad(\"/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/01.all.cluster/split_diff/pancancer.split.$pp.0831.h5ad\")
sc.set_figure_params(figsize=(15, 15))

";
open IN1,"/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/05.diff_six/diff/find_diff/merge.all.diffgene.xls";

$n=-1;
while(<IN1>){
chomp;
@a=split;if($a[-1] ne "diff.$tt"){next};
$n++;$ha{$a[1]}=1;
$m=int($n/50);
$hb{$m}.="\"$a[1]\"\,";
}

open OUT1,">/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/05.diff_six/diff/find_diff/find_plot/$tt/find.$tt.sh";

for $i(sort {$a<=>$b} keys %hb){
print  OUT0 "
marker_genes_dict=[$hb{$i}]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby=\"groups_pri\", var_group_labels=\"$tt\",
                  save=\"$tt.markergene.$i.rawpdf\")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby=\"groups_pri\", colorbar_title=\"mean z-score\",use_raw=False, vmin=-1, vmax=1, cmap=\"RdBu_r\", save=\"$tt.markergene.$i.xxx.pdf\")
"
}

print  OUT1 "/jdfssz1/ST_TSCBI/P22Z10200N0433/USER/zhangzhao/software/anaconda3/bin/python /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/05.diff_six/diff/find_diff/find_plot/$tt/find.$tt.py\n"
