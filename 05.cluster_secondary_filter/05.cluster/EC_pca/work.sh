perl -e 'print "cell\n";while(<>){chomp;@a=split(/\,/);if($a[-1] eq "7" || $a[-1] eq "11"){}else{print "$a[0]\n"}}' pancancer.ref.0723.final.EC.umap.obs.txt > EC.input &

perl -e '@x=(1000,2000);@y=(30);@z=(0.5,0.6,0.8,1,1.5);while(<>){chomp;$t=$_;$o=`cat cluster.py`;for($i=0;$i<@x;$i++){for($j=0;$j<@y;$j++){for($k=0;$k<@z;$k++){`mkdir -p $t/$x[$i]_$y[$j]_$z[$k]`;`ln -s /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/04.data_sub/01.secondary/pancancer.ref.0723.final.EC.h5ad $t/$x[$i]_$y[$j]_$z[$k]`}}}}' celltype.list 

perl -e '@x=(1000,2000);@y=(30);@z=(0.5,0.6,0.8,1,1.5);while(<>){chomp;$t=$_;for($i=0;$i<@x;$i++){for($j=0;$j<@y;$j++){for($k=0;$k<@z;$k++){$o=`cat cluster.py`;$p="$t/$x[$i]_$y[$j]_$z[$k]/cluster.$x[$i]_$y[$j]_$z[$k].py";open OUT,">$p";$o=~s/aaaa/$t/g;$o=~s/bbbb/$x[$i]/g;$o=~s/cccc/$y[$j]/g;$o=~s/dddd/$z[$k]/g;print OUT "$o";$p="$t/$x[$i]_$y[$j]_$z[$k]/cluster.$x[$i]_$y[$j]_$z[$k].sh";open OUT,">$p";print OUT "/jdfssz1/ST_TSCBI/P22Z10200N0433/USER/zhangzhao/software/anaconda3/bin/python /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/05.cluster_secondary_filter/05.cluster/EC_pca/$t/$x[$i]_$y[$j]_$z[$k]/cluster.$x[$i]_$y[$j]_$z[$k].py\n"}}}}' celltype.list

perl -e '@x=(1000,2000);@y=(30);@z=(0.5,0.6,0.8,1,1.5);while(<>){chomp;$t=$_;for($i=0;$i<@x;$i++){for($j=0;$j<@y;$j++){for($k=0;$k<@z;$k++){print "cd $t/$x[$i]_$y[$j]_$z[$k]\nqsub -cwd -l vf=50G,num_proc=1  -P P22Z10200N0433 -binding linear:1 -q st.q cluster.$x[$i]_$y[$j]_$z[$k].sh\ncd -\n"}}}}' celltype.list |sh

perl -e '@x=(1000,2000);@y=(30);@z=(0.5,0.6,0.8,1,1.5);while(<>){chomp;$t=$_;for($i=0;$i<@x;$i++){for($j=0;$j<@y;$j++){for($k=0;$k<@z;$k++){print "mkdir ../sub_EC_cluster/$t/$x[$i]_$y[$j]_$z[$k]\ncp $t/$x[$i]_$y[$j]_$z[$k]/figures/*png ../sub_EC_cluster/$t/$x[$i]_$y[$j]_$z[$k]\n"}}}}' celltype.list |sh


