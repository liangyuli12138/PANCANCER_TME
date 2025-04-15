perl -e '$t=<>;while(<>){chomp;@a=split(/,/);$ha{$a[-1]}.="$_\n"};for $i(keys %ha){open OUT,">split.$i.at";print OUT "$t","$ha{$i}"}' merge.all.all.level.filterB.sev.at
perl -e '<>;while(<>){chomp;@a=split(/,/);$ha{$a[-1]}.="$a[0]\n"};for $i(keys %ha){open OUT,">split.$i.input";print OUT "cell\n","$ha{$i}"}' merge.all.all.level.filterB.sev.at

perl -e '$t=`cat tmp.py`;while(<>){chomp;$b=$_;$p=$t;$p=~s/aaaa/$b/g;open OUT,">split.$b.py";print OUT "$p";open OUT,">split.$b.sh";print OUT "/jdfssz1/ST_TSCBI/P22Z10200N0433/USER/zhangzhao/software/anaconda3/bin/python split.$b.py\n";}' cell.list

perl -e 'while(<>){chomp;$b=$_;print "qsub -cwd -l vf=100G,num_proc=1  -P P22Z10200N0433 -binding linear:1 -q st.q split.$b.sh\n"}' cell.list |sh

