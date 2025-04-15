perl -e 'while(<>){chomp;@a=split(/,/);if($a[15] eq "3"){next};print "$a[0]\n"}' pancancer.ref.0807.final.CD4_pca.umap.obs.txt > filter.input

perl -e '@x=(2000);@y=(20);@z=(0.6,0.8,1,1.2,1.5,1.8,2,2.5,3,3.5);while(<>){chomp;$t=$_;$o=`cat cluster.py`;for($i=0;$i<@x;$i++){for($j=0;$j<@y;$j++){for($k=0;$k<@z;$k++){`mkdir -p $t/$x[$i]_$y[$j]_$z[$k]`;`ln -s /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/05.cluster_secondary_filter/04.data.filter_round3/CD4/pancancer.ref.0804.final.CD4.umap.h5ad $t/$x[$i]_$y[$j]_$z[$k]`}}}}' celltype.list

perl -e '@x=(2000);@y=(20);@z=(0.6,0.8,1,1.2,1.5,1.8,2,2.5,3,3.5);while(<>){chomp;$t=$_;for($i=0;$i<@x;$i++){for($j=0;$j<@y;$j++){for($k=0;$k<@z;$k++){$o=`cat cluster.py`;$p="$t/$x[$i]_$y[$j]_$z[$k]/cluster.$x[$i]_$y[$j]_$z[$k].py";open OUT,">$p";$o=~s/aaaa/$t/g;$o=~s/bbbb/$x[$i]/g;$o=~s/cccc/$y[$j]/g;$o=~s/dddd/$z[$k]/g;print OUT "$o";$p="$t/$x[$i]_$y[$j]_$z[$k]/cluster.$x[$i]_$y[$j]_$z[$k].sh";open OUT,">$p";print OUT "/jdfssz1/ST_TSCBI/P22Z10200N0433/USER/zhangzhao/software/anaconda3/bin/python /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/05.cluster_secondary_filter/05.cluster/CD4_pca_filter/$t/$x[$i]_$y[$j]_$z[$k]/cluster.$x[$i]_$y[$j]_$z[$k].py\n"}}}}' celltype.list

perl -e '@x=(2000);@y=(20);@z=(0.6,0.8,1,1.2,1.5,1.8,2,2.5,3,3.5);while(<>){chomp;$t=$_;for($i=0;$i<@x;$i++){for($j=0;$j<@y;$j++){for($k=0;$k<@z;$k++){print "cd $t/$x[$i]_$y[$j]_$z[$k]\nqsub -cwd -l vf=100G,num_proc=1  -P P22Z10200N0433 -binding linear:1 -q st.q cluster.$x[$i]_$y[$j]_$z[$k].sh\ncd -\n"}}}}' celltype.list |sh

for i in CD4_pca_filter.cluster/*;do for j in $i/*;do echo mv $j CD4_pca_filter.cluster/`basename $i`.`basename $j`;done;done|sh

