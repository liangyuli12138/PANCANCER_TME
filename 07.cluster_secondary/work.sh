perl -e 'while(<>){chomp;@a=split;print "mkdir $a[0]\n"}' sub.list |sh

perl -e '@x=(500,1000,1500,2000);@y=(10,15,20,25);@z=(0.5,1);while(<>){chomp;$t=$_;$o=`cat cluster.py`;$p=`cat marker_py/$t.py`;for($i=0;$i<@x;$i++){for($j=0;$j<@y;$j++){for($k=0;$k<@z;$k++){`mkdir -p $t/$x[$i]_$y[$j]_$z[$k]`;`ln -s /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/04.data_sub/01.secondary/pancancer.ref.0723.final.h5ad $t/$x[$i]_$y[$j]_$z[$k]`}}}}' sub.list 

perl -e '@x=(500,1000,1500,2000);@y=(10,15,20,25);@z=(0.5,1);while(<>){chomp;$t=$_;for($i=0;$i<@x;$i++){for($j=0;$j<@y;$j++){for($k=0;$k<@z;$k++){$o=`cat cluster.py marker_py/$t.py`;$p="$t/$x[$i]_$y[$j]_$z[$k]/cluster.$x[$i]_$y[$j]_$z[$k].py";open OUT,">$p";$o=~s/aaaa/$t/g;$o=~s/bbbb/$x[$i]/g;$o=~s/cccc/$y[$j]/g;$o=~s/dddd/$z[$k]/g;print OUT "$o";$p="$t/$x[$i]_$y[$j]_$z[$k]/cluster.$x[$i]_$y[$j]_$z[$k].sh";open OUT,">$p";print OUT "/jdfssz1/ST_TSCBI/P22Z10200N0433/USER/zhangzhao/software/anaconda3/bin/python /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/07.cluster_secondary/$t/$x[$i]_$y[$j]_$z[$k]/cluster.$x[$i]_$y[$j]_$z[$k].py\n"}}}}' sub.list

perl -e '@x=(500,1000,1500,2000);@y=(10,15,20,25);@z=(0.5,1);while(<>){chomp;$t=$_;for($i=0;$i<@x;$i++){for($j=0;$j<@y;$j++){for($k=0;$k<@z;$k++){print "cd $t/$x[$i]_$y[$j]_$z[$k]\nqsub -cwd -l vf=100G,num_proc=1  -P P22Z10200N0433 -binding linear:1 -q st.q cluster.$x[$i]_$y[$j]_$z[$k].sh\ncd -\n"}}}}' sub.list |sh

les sub.list |while read i;do echo -e "cp $i/*/*/*png all_cluster/$i";done|sh &

perl -e '@x=(1000,1500);@y=(10,15);@z=(0.5,0.8,1);@q=(0.1,0.3,0.6,1,2);while(<>){chomp;$t=$_;$o=`cat cluster.fib.py`;$p=`cat marker_py/$t.py`;for($i=0;$i<@x;$i++){for($j=0;$j<@y;$j++){for($k=0;$k<@z;$k++){for($l=0;$l<@q;$l++){`mkdir -p $t/$x[$i]_$y[$j]_$z[$k]_$q[$l]`;`ln -s /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/04.data_sub/01.secondary/pancancer.ref.0723.final.h5ad $t/$x[$i]_$y[$j]_$z[$k]_$q[$l]`}}}}}' sub2.list 

perl -e '@x=(1000,1500);@y=(10,15);@z=(0.5,0.8,1);@q=(0.1,0.3,0.6,1,2);while(<>){chomp;$t=$_;for($i=0;$i<@x;$i++){for($j=0;$j<@y;$j++){for($k=0;$k<@z;$k++){for($l=0;$l<@q;$l++){$o=`cat cluster.fib.py marker_py/$t.py`;$p="$t/$x[$i]_$y[$j]_$z[$k]_$q[$l]/cluster.$x[$i]_$y[$j]_$z[$k]_$q[$l].py";open OUT,">$p";$o=~s/aaaa/$t/g;$o=~s/bbbb/$x[$i]/g;$o=~s/cccc/$y[$j]/g;$o=~s/dddd/$z[$k]/g;$o=~s/eeee/$q[$l]/g;print OUT "$o";$p="$t/$x[$i]_$y[$j]_$z[$k]_$q[$l]/cluster.$x[$i]_$y[$j]_$z[$k]_$q[$l].sh";open OUT,">$p";print OUT "/jdfssz1/ST_TSCBI/P22Z10200N0433/USER/zhangzhao/software/anaconda3/bin/python /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/07.cluster_secondary/$t/$x[$i]_$y[$j]_$z[$k]_$q[$l]/cluster.$x[$i]_$y[$j]_$z[$k]_$q[$l].py\n"}}}}}' sub2.list 

perl -e '@x=(1000,1500);@y=(10,15);@z=(0.5,0.8,1);@q=(0.1,0.3,0.6,1,2);while(<>){chomp;$t=$_;for($i=0;$i<@x;$i++){for($j=0;$j<@y;$j++){for($k=0;$k<@z;$k++){for($l=0;$l<@q;$l++){print "cd $t/$x[$i]_$y[$j]_$z[$k]_$q[$l]\nqsub -cwd -l vf=50G,num_proc=1  -P P22Z10200N0433 -binding linear:1 -q st.q cluster.$x[$i]_$y[$j]_$z[$k]_$q[$l].sh\ncd -\n"}}}}} ' sub2.list |sh

