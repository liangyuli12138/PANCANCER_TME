perl -e 'while(<>){chomp;@a=split(/\//);print "mkdir -p diff/$a[9]/$a[10]\n"}'  all.obs.list |sh

perl -e 'while(<>){chomp;open IN,$_;$n=0;$t=$_;@z=split(/\//,$t);undef %ha;while(<IN>){chomp;@a=split(/\,/);if(!exists $ha{$a[-1]}){$n++;$ha{$a[-1]}=1}};$n=$n-1;$x=$t;$t=~s/obs\.txt/h5ad/;for($i=0;$i<$n;$i++){print "/jdfssz1/ST_TSCBI/P22Z10200N0433/USER/zhangzhao/software/anaconda3/bin/python /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/05.cluster_secondary_filter/03.diff/findmarker.py $t $i /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/07.cluster_secondary/all_diff/diff/$z[9]/$z[10]/diff.$i.csv\n"}}'  all.obs.list > all.diff.sh &

perl -e 'while(<>){chomp;@a=split(/\//);print "mkdir -p all_sub_diff/$a[9]/$a[10]\n"}'  all.obs.list |sh

perl -e 'while(<>){chomp;$a=$_;$b=`head -n 501 $a`;$a=~s/diff/all_sub_diff/;open OUT,">$a";print OUT "$b"}' all.diff.list &

perl -e 'while(<>){chomp;@a=split(/\//);print "perl /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/05.cluster_secondary/all_diff/overlap.pl /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancancer/08.filter_by_gene/13.ref_20230604/08.all.h5ad_ref/anndata_0.8/pancancer.ref.0621.raw.obs.csv $_ > /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/07.cluster_secondary/all_diff/all_sub_diff/$a[9]/$a[10]/overlap.stat.ori.txt\n"}' all.obs.list |sh &

perl -e 'while(<>){chomp;@a=split(/\//);print "perl overlap.pl $_ $_ > /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/07.cluster_secondary/all_diff/all_sub_diff/$a[9]/$a[10]/overlap.stat.new.txt\n"}' all.obs.list |sh &

