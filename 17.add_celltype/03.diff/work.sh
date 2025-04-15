perl -e 'while(<>){chomp;@a=split(/\//);print "mkdir -p $a[10]/$a[11]\n"}'  all.obs.list |sh

perl -e 'while(<>){chomp;open IN,$_;$n=0;$t=$_;@z=split(/\//,$t);undef %ha;while(<IN>){chomp;@a=split(/\,/);if(!exists $ha{$a[-1]}){$n++;$ha{$a[-1]}=1}};$n=$n-1;$x=$t;$t=~s/obs\.txt/h5ad/;for($i=0;$i<$n;$i++){print "python /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/05.cluster_secondary/all_diff/findmarker.py $t $i /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/17.add_celltype/03.diff/$z[10]/$z[11]/diff.$i.csv\n"}}'  all.obs.list > all.diff.sh &

perl -e 'while(<>){chomp;@a=split(/\//);print "perl /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/05.cluster_secondary/all_diff/overlap.pl /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancancer/08.filter_by_gene/13.ref_20230604/08.all.h5ad_ref/anndata_0.8/pancancer.ref.0621.raw.obs.csv $_ > /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/05.cluster_secondary/all_diff/$a[9]/$a[10]/overlap.stat.txt\n"}' all.obs.list > all.overlap.sh

