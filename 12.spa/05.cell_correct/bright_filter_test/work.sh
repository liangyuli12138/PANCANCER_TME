perl -e 'while(<>){chomp;$a=$_;$o=`cat test.py`;$o=~s/aaaa/$a/g;open OUT,">shell/bright.$a.py";print OUT "$o";open OUT,">shell/bright.$a.sh";print OUT "/ldfssz1/ST_OCEAN/USER/liaoshangfeng/software/anaconda3/envs/cell2loc_env/bin/python /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/12.spa/05.cell_correct/bright_filter_test/shell/bright.$a.py\n"}' sample_list 

