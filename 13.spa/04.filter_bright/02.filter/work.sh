perl -e 'while(<>){chomp;@a=split;$o=`cat test.filter.py`;$o=~s/aaaa/$a[0]/g;$o=~s/bbbb/$a[1]/g;open OUT,">shell/filter.$a[0].py";print OUT "$o";open OUT,">shell/filter.$a[0].sh";print OUT "/ldfssz1/ST_OCEAN/USER/liaoshangfeng/software/anaconda3/envs/cell2loc_env/bin/python filter.$a[0].py\n"}' filter.parameter 

