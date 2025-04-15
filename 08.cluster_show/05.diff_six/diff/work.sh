ls  /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/01.all.cluster/split_diff/*obs* > obs.list

perl -e 'while(<>){chomp;open IN,$_;$n=0;$t=$_;@z=split(/\//,$t);undef %ha;<IN>;while(<IN>){chomp;@a=split(/\,/);if(!exists $ha{$a[-5]}){$n++;$ha{$a[-5]}=1}};$n=$n-1;$x=$t;$t=~s/obs\.csv/h5ad/;for $i(keys %ha){print "/jdfssz1/ST_TSCBI/P22Z10200N0433/USER/zhangzhao/software/anaconda3/bin/python /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/05.diff_six/diff/findmarker.py $t $i /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/05.diff_six/diff/diff/diff.$i.csv\n"}}' obs.list > all.diff.sh

