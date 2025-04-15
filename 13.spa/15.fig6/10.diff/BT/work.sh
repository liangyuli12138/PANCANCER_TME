perl -e 'while(<>){chomp;@a=split(/,/);print "$a[-1]\n"}' pancancer.icar.bt.bg.cell.obs > list

perl -e 'while(<>){chomp;@a=split;$ha{$a[0]}=1};for $i(keys %ha){for $j(keys %ha){if($i eq $j){next};@x=split(/\./,$i);@y=split(/\./,$j);if($x[1] ne $y[1]){next};print "/hwfssz4/BC_PUB/Software/07.User-defined/03.Animal_Plant/wubin/mamba/bin/python /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/15.fig6/10.diff/BT/findmarker.py /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/15.fig6/10.diff/pancancer.icar.bt.bg.cell.h5ad $i $j /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/15.fig6/10.diff/BT/diff/diff.$i.$j.csv\n"}}' list > all.sh

