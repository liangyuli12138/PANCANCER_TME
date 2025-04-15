perl -e 'while(<>){chomp;@a=split;$o=`cat CellBin_RCTD_plot.hg.all.R`;$o=~s/bbbb/$a[1].$a[0]/g;$o=~s/aaaa/$a[0]/g;open OUT,">shell/plot.$a[0].R";print OUT "$o";open OUT,">shell/plot.$a[0].sh";print OUT "/ldfssz1/ST_OCEAN/USER/liaoshangfeng/software/anaconda3/envs/R411/bin/Rscript plot.$a[0].R\n";}' all.sn.list

